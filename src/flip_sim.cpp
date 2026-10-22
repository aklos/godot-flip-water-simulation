#include "flip_sim.h"
#include "flip_effectors.h"

#include <godot_cpp/classes/engine.hpp>
#include <godot_cpp/classes/os.hpp>
#include <godot_cpp/classes/physics_direct_space_state2d.hpp>
#include <godot_cpp/classes/physics_server2d.hpp>
#include <godot_cpp/classes/rendering_server.hpp>
#include <godot_cpp/classes/scene_tree.hpp>
#include <godot_cpp/classes/shader.hpp>
#include <godot_cpp/classes/sprite2d.hpp>
#include <godot_cpp/classes/sub_viewport.hpp>
#include <godot_cpp/classes/viewport_texture.hpp>
#include <godot_cpp/classes/worker_thread_pool.hpp>
#include <godot_cpp/classes/world2d.hpp>
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/error_macros.hpp>
#include <godot_cpp/core/object.hpp>
#include <godot_cpp/variant/callable_method_pointer.hpp>
#include <godot_cpp/variant/dictionary.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace godot;

static const double DENSITY = 1000.0;
static const double OVER_RELAXATION = 1.9;
// How fast (in cells per second, per unit of excess density) over-compressed regions are pushed apart.
static const double DRIFT_COMPENSATION = 5.0;
static const int MAX_QUERY_RESULTS = 1024;
static const int RENDER_FLOATS_PER_INSTANCE = 12; // 2D transform (8) + color (4).

static const char *PARTICLE_SHADER = R"(
shader_type canvas_item;

void fragment() {
	float d = length(UV - vec2(0.5)) * 2.0;
	COLOR.a *= 1.0 - smoothstep(0.85, 1.0, d);
}
)";

// Additively splats a soft kernel per particle: R accumulates density, G accumulates foam.
static const char *SPLAT_SHADER = R"(
shader_type canvas_item;
render_mode blend_add, unshaded;

void fragment() {
	float d = length(UV - vec2(0.5)) * 2.0;
	float w = clamp(1.0 - d * d, 0.0, 1.0);
	w *= w * 0.5;
	COLOR = vec4(w, w * COLOR.r, 0.0, 1.0);
}
)";

// Thresholds the splatted density into a smooth liquid surface.
static const char *SURFACE_SHADER = R"(
shader_type canvas_item;

uniform vec4 water_color : source_color = vec4(0.2, 0.5, 1.0, 0.9);
uniform vec4 foam_color : source_color = vec4(0.85, 0.95, 1.0, 1.0);
uniform float threshold = 0.25;

void fragment() {
	vec4 t = texture(TEXTURE, UV);
	float a = smoothstep(threshold, threshold + 0.08, t.r);
	float foam = clamp(t.g / max(t.r, 0.0001), 0.0, 1.0);
	vec4 c = mix(water_color, foam_color, foam);
	COLOR = vec4(c.rgb, c.a * a);
}
)";

template <typename T>
static void shift_grid(std::vector<T> &a, int nx, int ny, int dx, int dy, T fill) {
	std::vector<T> out(a.size(), fill);
	for (int j = 0; j < ny; j++) {
		int sj = j + dy;
		if (sj < 0 || sj >= ny) {
			continue;
		}
		for (int i = 0; i < nx; i++) {
			int si = i + dx;
			if (si >= 0 && si < nx) {
				out[i + (nx * j)] = a[si + (nx * sj)];
			}
		}
	}
	a.swap(out);
}

static double segment_distance(const Vector2 &p, const Vector2 &a, const Vector2 &b) {
	Vector2 ab = b - a;
	double len2 = ab.length_squared();
	double t = len2 > 0.0 ? std::clamp(static_cast<double>((p - a).dot(ab)) / len2, 0.0, 1.0) : 0.0;
	return (p - (a + ab * t)).length();
}

static Ref<ShaderMaterial> make_material(const char *code) {
	Ref<Shader> shader;
	shader.instantiate();
	shader->set_code(code);
	Ref<ShaderMaterial> material;
	material.instantiate();
	material->set_shader(shader);
	return material;
}

void FLIPWaterSimulation::_bind_methods() {
	ClassDB::bind_method(D_METHOD("set_grid_size", "size"), &FLIPWaterSimulation::set_grid_size);
	ClassDB::bind_method(D_METHOD("get_grid_size"), &FLIPWaterSimulation::get_grid_size);
	ClassDB::bind_method(D_METHOD("set_cell_size", "size"), &FLIPWaterSimulation::set_cell_size);
	ClassDB::bind_method(D_METHOD("get_cell_size"), &FLIPWaterSimulation::get_cell_size);
	ClassDB::bind_method(D_METHOD("set_gravity", "gravity"), &FLIPWaterSimulation::set_gravity);
	ClassDB::bind_method(D_METHOD("get_gravity"), &FLIPWaterSimulation::get_gravity);
	ClassDB::bind_method(D_METHOD("set_flip_ratio", "flip_ratio"), &FLIPWaterSimulation::set_flip_ratio);
	ClassDB::bind_method(D_METHOD("get_flip_ratio"), &FLIPWaterSimulation::get_flip_ratio);
	ClassDB::bind_method(D_METHOD("set_max_speed", "max_speed"), &FLIPWaterSimulation::set_max_speed);
	ClassDB::bind_method(D_METHOD("get_max_speed"), &FLIPWaterSimulation::get_max_speed);
	ClassDB::bind_method(D_METHOD("set_solver_iterations", "iterations"), &FLIPWaterSimulation::set_solver_iterations);
	ClassDB::bind_method(D_METHOD("get_solver_iterations"), &FLIPWaterSimulation::get_solver_iterations);
	ClassDB::bind_method(D_METHOD("set_substeps", "substeps"), &FLIPWaterSimulation::set_substeps);
	ClassDB::bind_method(D_METHOD("get_substeps"), &FLIPWaterSimulation::get_substeps);
	ClassDB::bind_method(D_METHOD("set_max_particles", "max_particles"), &FLIPWaterSimulation::set_max_particles);
	ClassDB::bind_method(D_METHOD("get_max_particles"), &FLIPWaterSimulation::get_max_particles);
	ClassDB::bind_method(D_METHOD("set_edge_mode", "side", "mode"), &FLIPWaterSimulation::set_edge_mode);
	ClassDB::bind_method(D_METHOD("get_edge_mode", "side"), &FLIPWaterSimulation::get_edge_mode);
	ClassDB::bind_method(D_METHOD("set_auto_update", "enabled"), &FLIPWaterSimulation::set_auto_update);
	ClassDB::bind_method(D_METHOD("is_auto_update"), &FLIPWaterSimulation::is_auto_update);
	ClassDB::bind_method(D_METHOD("set_use_threads", "enabled"), &FLIPWaterSimulation::set_use_threads);
	ClassDB::bind_method(D_METHOD("is_using_threads"), &FLIPWaterSimulation::is_using_threads);
	ClassDB::bind_method(D_METHOD("set_follow_target", "path"), &FLIPWaterSimulation::set_follow_target);
	ClassDB::bind_method(D_METHOD("get_follow_target"), &FLIPWaterSimulation::get_follow_target);
	ClassDB::bind_method(D_METHOD("set_collision_mask", "mask"), &FLIPWaterSimulation::set_collision_mask);
	ClassDB::bind_method(D_METHOD("get_collision_mask"), &FLIPWaterSimulation::get_collision_mask);
	ClassDB::bind_method(D_METHOD("set_debug_draw", "enabled"), &FLIPWaterSimulation::set_debug_draw);
	ClassDB::bind_method(D_METHOD("is_debug_draw"), &FLIPWaterSimulation::is_debug_draw);

	ClassDB::bind_method(D_METHOD("set_render_mode", "mode"), &FLIPWaterSimulation::set_render_mode);
	ClassDB::bind_method(D_METHOD("get_render_mode"), &FLIPWaterSimulation::get_render_mode);
	ClassDB::bind_method(D_METHOD("set_water_color", "color"), &FLIPWaterSimulation::set_water_color);
	ClassDB::bind_method(D_METHOD("get_water_color"), &FLIPWaterSimulation::get_water_color);
	ClassDB::bind_method(D_METHOD("set_foam_color", "color"), &FLIPWaterSimulation::set_foam_color);
	ClassDB::bind_method(D_METHOD("get_foam_color"), &FLIPWaterSimulation::get_foam_color);
	ClassDB::bind_method(D_METHOD("set_foam_speed", "speed"), &FLIPWaterSimulation::set_foam_speed);
	ClassDB::bind_method(D_METHOD("get_foam_speed"), &FLIPWaterSimulation::get_foam_speed);
	ClassDB::bind_method(D_METHOD("set_particle_scale", "scale"), &FLIPWaterSimulation::set_particle_scale);
	ClassDB::bind_method(D_METHOD("get_particle_scale"), &FLIPWaterSimulation::get_particle_scale);
	ClassDB::bind_method(D_METHOD("set_surface_resolution", "resolution"), &FLIPWaterSimulation::set_surface_resolution);
	ClassDB::bind_method(D_METHOD("get_surface_resolution"), &FLIPWaterSimulation::get_surface_resolution);

	ClassDB::bind_method(D_METHOD("get_particle_radius"), &FLIPWaterSimulation::get_particle_radius);
	ClassDB::bind_method(D_METHOD("get_domain_rect"), &FLIPWaterSimulation::get_domain_rect);
	ClassDB::bind_method(D_METHOD("local_to_cell", "position"), &FLIPWaterSimulation::local_to_cell);
	ClassDB::bind_method(D_METHOD("cell_to_local", "cell"), &FLIPWaterSimulation::cell_to_local);
	ClassDB::bind_method(D_METHOD("set_cell", "cell", "type"), &FLIPWaterSimulation::set_cell);
	ClassDB::bind_method(D_METHOD("get_cell", "cell"), &FLIPWaterSimulation::get_cell);
	ClassDB::bind_method(D_METHOD("fill_cells", "rect", "type"), &FLIPWaterSimulation::fill_cells);
	ClassDB::bind_method(D_METHOD("clear"), &FLIPWaterSimulation::clear);
	ClassDB::bind_method(D_METHOD("shift", "cells"), &FLIPWaterSimulation::shift);
	ClassDB::bind_method(D_METHOD("refresh_colliders"), &FLIPWaterSimulation::refresh_colliders);
	ClassDB::bind_method(D_METHOD("update", "delta"), &FLIPWaterSimulation::update);

	ClassDB::bind_method(D_METHOD("get_particle_count"), &FLIPWaterSimulation::get_particle_count);
	ClassDB::bind_method(D_METHOD("get_particle_positions"), &FLIPWaterSimulation::get_particle_positions);
	ClassDB::bind_method(D_METHOD("get_particle_velocities"), &FLIPWaterSimulation::get_particle_velocities);
	ClassDB::bind_method(D_METHOD("get_cell_types"), &FLIPWaterSimulation::get_cell_types);
	ClassDB::bind_method(D_METHOD("get_pressure"), &FLIPWaterSimulation::get_pressure);
	ClassDB::bind_method(D_METHOD("add_particle", "position", "velocity"), &FLIPWaterSimulation::add_particle, DEFVAL(Vector2()));
	ClassDB::bind_method(D_METHOD("apply_impulse", "center", "radius", "velocity"), &FLIPWaterSimulation::apply_impulse);
	ClassDB::bind_method(D_METHOD("remove_particles_in_circle", "center", "radius"), &FLIPWaterSimulation::remove_particles_in_circle);
	ClassDB::bind_method(D_METHOD("count_particles_in_circle", "center", "radius"), &FLIPWaterSimulation::count_particles_in_circle);
	ClassDB::bind_method(D_METHOD("count_particles_in_rect", "rect"), &FLIPWaterSimulation::count_particles_in_rect);

	ADD_GROUP("Grid", "");
	ADD_PROPERTY(PropertyInfo(Variant::VECTOR2I, "grid_size", PROPERTY_HINT_NONE, "suffix:cells"), "set_grid_size", "get_grid_size");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "cell_size", PROPERTY_HINT_RANGE, "1,256,0.1,or_greater,suffix:px"), "set_cell_size", "get_cell_size");
	ADD_PROPERTYI(PropertyInfo(Variant::INT, "edge_left", PROPERTY_HINT_ENUM, "Wall,Open"), "set_edge_mode", "get_edge_mode", SIDE_LEFT);
	ADD_PROPERTYI(PropertyInfo(Variant::INT, "edge_top", PROPERTY_HINT_ENUM, "Wall,Open"), "set_edge_mode", "get_edge_mode", SIDE_TOP);
	ADD_PROPERTYI(PropertyInfo(Variant::INT, "edge_right", PROPERTY_HINT_ENUM, "Wall,Open"), "set_edge_mode", "get_edge_mode", SIDE_RIGHT);
	ADD_PROPERTYI(PropertyInfo(Variant::INT, "edge_bottom", PROPERTY_HINT_ENUM, "Wall,Open"), "set_edge_mode", "get_edge_mode", SIDE_BOTTOM);

	ADD_GROUP("Simulation", "");
	ADD_PROPERTY(PropertyInfo(Variant::VECTOR2, "gravity", PROPERTY_HINT_NONE, "suffix:px/s²"), "set_gravity", "get_gravity");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "flip_ratio", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_flip_ratio", "get_flip_ratio");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "max_speed", PROPERTY_HINT_RANGE, "0,5000,1,or_greater,suffix:px/s"), "set_max_speed", "get_max_speed");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "solver_iterations", PROPERTY_HINT_RANGE, "1,200,1"), "set_solver_iterations", "get_solver_iterations");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "substeps", PROPERTY_HINT_RANGE, "1,16,1"), "set_substeps", "get_substeps");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "max_particles", PROPERTY_HINT_RANGE, "0,200000,1,or_greater"), "set_max_particles", "get_max_particles");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "auto_update"), "set_auto_update", "is_auto_update");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "use_threads"), "set_use_threads", "is_using_threads");

	ADD_GROUP("Scrolling", "");
	ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH, "follow_target", PROPERTY_HINT_NODE_PATH_VALID_TYPES, "Node2D"), "set_follow_target", "get_follow_target");

	ADD_GROUP("Collision", "");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "collision_mask", PROPERTY_HINT_LAYERS_2D_PHYSICS), "set_collision_mask", "get_collision_mask");

	ADD_GROUP("Rendering", "");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "render_mode", PROPERTY_HINT_ENUM, "None,Particles,Surface"), "set_render_mode", "get_render_mode");
	ADD_PROPERTY(PropertyInfo(Variant::COLOR, "water_color"), "set_water_color", "get_water_color");
	ADD_PROPERTY(PropertyInfo(Variant::COLOR, "foam_color"), "set_foam_color", "get_foam_color");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "foam_speed", PROPERTY_HINT_RANGE, "0,5000,1,or_greater,suffix:px/s"), "set_foam_speed", "get_foam_speed");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "particle_scale", PROPERTY_HINT_RANGE, "0.1,4,0.01"), "set_particle_scale", "get_particle_scale");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "surface_resolution", PROPERTY_HINT_RANGE, "0.1,1,0.05"), "set_surface_resolution", "get_surface_resolution");

	ADD_GROUP("Debug", "");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "debug_draw"), "set_debug_draw", "is_debug_draw");

	ADD_SIGNAL(MethodInfo("region_exposed", PropertyInfo(Variant::RECT2I, "rect")));

	BIND_ENUM_CONSTANT(AIR_CELL);
	BIND_ENUM_CONSTANT(SOLID_CELL);
	BIND_ENUM_CONSTANT(FLUID_CELL);
	BIND_ENUM_CONSTANT(EDGE_WALL);
	BIND_ENUM_CONSTANT(EDGE_OPEN);
	BIND_ENUM_CONSTANT(RENDER_NONE);
	BIND_ENUM_CONSTANT(RENDER_PARTICLES);
	BIND_ENUM_CONSTANT(RENDER_SURFACE);
}

FLIPWaterSimulation::FLIPWaterSimulation() {
	grid_size = Vector2i(64, 36);
	cell_size = 20.0;
	gravity = Vector2(0.0, 980.0);
	flip_ratio = 0.9;
	max_speed = 0.0;
	solver_iterations = 30;
	substeps = 1;
	max_particles = 20000;
	for (int i = 0; i < 4; i++) {
		edge_modes[i] = EDGE_WALL;
	}
	auto_update = true;
	use_threads = true;
	collision_mask = 1;
	debug_draw = false;

	render_mode = RENDER_PARTICLES;
	water_color = Color(0.2, 0.5, 1.0, 0.9);
	foam_color = Color(0.85, 0.95, 1.0, 1.0);
	foam_speed = 600.0;
	particle_scale = 1.0;
	surface_resolution = 0.5;

	rng.seed(12345);
	static_colliders_dirty = true;
	job_fn = nullptr;
	job_count = 0;
	job_chunks = 0;
	solve_color = 0;
	solve_cp = 0.0;
	surface_viewport = nullptr;
	surface_root = nullptr;
	surface_sprite = nullptr;
	render_capacity = 0;

	_rebuild_grid();
}

FLIPWaterSimulation::~FLIPWaterSimulation() {
	_free_render_resources();
	if (query_shape.is_valid()) {
		PhysicsServer2D::get_singleton()->free_rid(query_shape);
	}
}

void FLIPWaterSimulation::_notification(int p_what) {
	bool editor = Engine::get_singleton()->is_editor_hint();
	switch (p_what) {
		case NOTIFICATION_ENTER_TREE:
			static_colliders_dirty = true;
			break;
		case NOTIFICATION_READY:
			if (!editor) {
				_ensure_render_resources();
				_update_render();
				set_physics_process(true);
			}
			break;
		case NOTIFICATION_PHYSICS_PROCESS:
			if (auto_update && !editor) {
				_follow();
				update(get_physics_process_delta_time());
			}
			break;
		case NOTIFICATION_DRAW: {
			if (!editor && !debug_draw) {
				break;
			}
			Rect2 domain = get_domain_rect();
			Color wall_color(0.5, 0.5, 0.55, 0.35);
			double h = cell_size;
			if (edge_modes[SIDE_LEFT] == EDGE_WALL) {
				draw_rect(Rect2(0, 0, h, domain.size.y), wall_color);
			}
			if (edge_modes[SIDE_RIGHT] == EDGE_WALL) {
				draw_rect(Rect2(domain.size.x - h, 0, h, domain.size.y), wall_color);
			}
			if (edge_modes[SIDE_TOP] == EDGE_WALL) {
				draw_rect(Rect2(0, 0, domain.size.x, h), wall_color);
			}
			if (edge_modes[SIDE_BOTTOM] == EDGE_WALL) {
				draw_rect(Rect2(0, domain.size.y - h, domain.size.x, h), wall_color);
			}
			if (!editor) {
				for (int j = 0; j < f_num_y; j++) {
					for (int i = 0; i < f_num_x; i++) {
						int idx = i + (f_num_x * j);
						Rect2 r(i * h, j * h, h, h);
						if (dynamic_solid[idx]) {
							draw_rect(r, Color(1.0, 0.5, 0.2, 0.6));
						} else if (user_solid[idx] || static_solid[idx]) {
							draw_rect(r, Color(0.5, 0.5, 0.55, 0.6));
						} else if (cell_types[idx] == FLUID_CELL) {
							draw_rect(r, Color(0.2, 0.5, 1.0, 0.15));
						}
					}
				}
			}
			draw_rect(domain, Color(0.35, 0.65, 1.0, 0.8), false, 2.0);
		} break;
	}
}

// Grid setup

void FLIPWaterSimulation::_rebuild_grid() {
	f_num_x = grid_size.x;
	f_num_y = grid_size.y;
	f_num_cells = f_num_x * f_num_y;
	f_inv_spacing = 1.0 / cell_size;
	p_radius = 0.3 * cell_size;
	p_rest_density = 0.0;

	cell_types.assign(f_num_cells, AIR_CELL);
	user_solid.assign(f_num_cells, 0);
	static_solid.assign(f_num_cells, 0);
	dynamic_solid.assign(f_num_cells, 0);
	solid_velocity.assign(f_num_cells, Vector2());
	u.assign(f_num_cells, 0.0);
	v.assign(f_num_cells, 0.0);
	du.assign(f_num_cells, 0.0);
	dv.assign(f_num_cells, 0.0);
	prev_u.assign(f_num_cells, 0.0);
	prev_v.assign(f_num_cells, 0.0);
	pressure.assign(f_num_cells, 0.0);
	p_densities.assign(f_num_cells, 0.0);
	s.assign(f_num_cells, 1.0);

	p_inv_spacing = 1.0 / (2.2 * p_radius);
	p_num_x = static_cast<int>(std::floor(f_num_x * cell_size * p_inv_spacing)) + 1;
	p_num_y = static_cast<int>(std::floor(f_num_y * cell_size * p_inv_spacing)) + 1;
	int p_num_cells = p_num_x * p_num_y;
	num_cell_particles.assign(p_num_cells, 0);
	first_cell_particle.assign(p_num_cells + 1, 0);

	p_positions.clear();
	p_velocities.clear();
	dynamic_cells.clear();
	static_colliders_dirty = true;
	_update_solids();
	queue_redraw();
}

bool FLIPWaterSimulation::_is_cell_in_bounds(int i, int j) const {
	return i >= 0 && i < f_num_x && j >= 0 && j < f_num_y;
}

bool FLIPWaterSimulation::_is_edge_wall(int i, int j) const {
	return (i == 0 && edge_modes[SIDE_LEFT] == EDGE_WALL) ||
			(j == 0 && edge_modes[SIDE_TOP] == EDGE_WALL) ||
			(i == f_num_x - 1 && edge_modes[SIDE_RIGHT] == EDGE_WALL) ||
			(j == f_num_y - 1 && edge_modes[SIDE_BOTTOM] == EDGE_WALL);
}

void FLIPWaterSimulation::_refresh_cell(int i, int j) {
	int idx = i + (f_num_x * j);
	bool solid = user_solid[idx] || static_solid[idx] || dynamic_solid[idx] || _is_edge_wall(i, j);
	s[idx] = solid ? 0.0 : 1.0;
	if (solid) {
		cell_types[idx] = SOLID_CELL;
	} else if (cell_types[idx] == SOLID_CELL) {
		cell_types[idx] = AIR_CELL;
	}
}

void FLIPWaterSimulation::_update_solids() {
	for (int j = 0; j < f_num_y; j++) {
		for (int i = 0; i < f_num_x; i++) {
			_refresh_cell(i, j);
		}
	}
}

// Properties

void FLIPWaterSimulation::set_grid_size(const Vector2i &p_size) {
	grid_size = Vector2i(std::max(p_size.x, 3), std::max(p_size.y, 3));
	_rebuild_grid();
}

Vector2i FLIPWaterSimulation::get_grid_size() const {
	return grid_size;
}

void FLIPWaterSimulation::set_cell_size(double p_size) {
	cell_size = std::max(p_size, 1.0);
	_rebuild_grid();
}

double FLIPWaterSimulation::get_cell_size() const {
	return cell_size;
}

void FLIPWaterSimulation::set_gravity(const Vector2 &p_gravity) {
	gravity = p_gravity;
}

Vector2 FLIPWaterSimulation::get_gravity() const {
	return gravity;
}

void FLIPWaterSimulation::set_flip_ratio(double p_flip_ratio) {
	flip_ratio = std::clamp(p_flip_ratio, 0.0, 1.0);
}

double FLIPWaterSimulation::get_flip_ratio() const {
	return flip_ratio;
}

void FLIPWaterSimulation::set_max_speed(double p_max_speed) {
	max_speed = std::max(p_max_speed, 0.0);
}

double FLIPWaterSimulation::get_max_speed() const {
	return max_speed;
}

void FLIPWaterSimulation::set_solver_iterations(int p_iterations) {
	solver_iterations = std::max(p_iterations, 1);
}

int FLIPWaterSimulation::get_solver_iterations() const {
	return solver_iterations;
}

void FLIPWaterSimulation::set_substeps(int p_substeps) {
	substeps = std::max(p_substeps, 1);
}

int FLIPWaterSimulation::get_substeps() const {
	return substeps;
}

void FLIPWaterSimulation::set_max_particles(int p_max) {
	max_particles = std::max(p_max, 0);
}

int FLIPWaterSimulation::get_max_particles() const {
	return max_particles;
}

void FLIPWaterSimulation::set_edge_mode(int p_side, EdgeMode p_mode) {
	ERR_FAIL_INDEX(p_side, 4);
	edge_modes[p_side] = p_mode;
	_update_solids();
	queue_redraw();
}

FLIPWaterSimulation::EdgeMode FLIPWaterSimulation::get_edge_mode(int p_side) const {
	ERR_FAIL_INDEX_V(p_side, 4, EDGE_WALL);
	return edge_modes[p_side];
}

void FLIPWaterSimulation::set_auto_update(bool p_enabled) {
	auto_update = p_enabled;
}

bool FLIPWaterSimulation::is_auto_update() const {
	return auto_update;
}

void FLIPWaterSimulation::set_use_threads(bool p_enabled) {
	use_threads = p_enabled;
}

bool FLIPWaterSimulation::is_using_threads() const {
	return use_threads;
}

void FLIPWaterSimulation::set_follow_target(const NodePath &p_path) {
	follow_target = p_path;
}

NodePath FLIPWaterSimulation::get_follow_target() const {
	return follow_target;
}

void FLIPWaterSimulation::set_collision_mask(uint32_t p_mask) {
	collision_mask = p_mask;
	static_colliders_dirty = true;
}

uint32_t FLIPWaterSimulation::get_collision_mask() const {
	return collision_mask;
}

void FLIPWaterSimulation::set_debug_draw(bool p_enabled) {
	debug_draw = p_enabled;
	queue_redraw();
}

bool FLIPWaterSimulation::is_debug_draw() const {
	return debug_draw;
}

void FLIPWaterSimulation::set_render_mode(RenderMode p_mode) {
	render_mode = p_mode;
	_update_render();
}

FLIPWaterSimulation::RenderMode FLIPWaterSimulation::get_render_mode() const {
	return render_mode;
}

void FLIPWaterSimulation::set_water_color(const Color &p_color) {
	water_color = p_color;
}

Color FLIPWaterSimulation::get_water_color() const {
	return water_color;
}

void FLIPWaterSimulation::set_foam_color(const Color &p_color) {
	foam_color = p_color;
}

Color FLIPWaterSimulation::get_foam_color() const {
	return foam_color;
}

void FLIPWaterSimulation::set_foam_speed(double p_speed) {
	foam_speed = std::max(p_speed, 0.0);
}

double FLIPWaterSimulation::get_foam_speed() const {
	return foam_speed;
}

void FLIPWaterSimulation::set_particle_scale(double p_scale) {
	particle_scale = std::max(p_scale, 0.01);
}

double FLIPWaterSimulation::get_particle_scale() const {
	return particle_scale;
}

void FLIPWaterSimulation::set_surface_resolution(double p_resolution) {
	surface_resolution = std::clamp(p_resolution, 0.05, 1.0);
}

double FLIPWaterSimulation::get_surface_resolution() const {
	return surface_resolution;
}

// Cells

double FLIPWaterSimulation::get_particle_radius() const {
	return p_radius;
}

Rect2 FLIPWaterSimulation::get_domain_rect() const {
	return Rect2(0, 0, f_num_x * cell_size, f_num_y * cell_size);
}

Vector2i FLIPWaterSimulation::local_to_cell(const Vector2 &p_pos) const {
	return Vector2i(static_cast<int>(std::floor(p_pos.x / cell_size)), static_cast<int>(std::floor(p_pos.y / cell_size)));
}

Vector2 FLIPWaterSimulation::cell_to_local(const Vector2i &p_cell) const {
	return Vector2(p_cell) * cell_size;
}

void FLIPWaterSimulation::_spawn_cell_particles(double x, double y) {
	for (int n = 0; n < 2; n++) {
		double offset = (cell_size / 3) + (n * (cell_size / 2));
		_add_particle_internal(Vector2(x + offset, y + offset), Vector2());
	}
}

void FLIPWaterSimulation::set_cell(const Vector2i &p_cell, CellType p_type) {
	if (!_is_cell_in_bounds(p_cell.x, p_cell.y)) {
		return;
	}
	int idx = p_cell.x + (f_num_x * p_cell.y);
	user_solid[idx] = p_type == SOLID_CELL;
	_refresh_cell(p_cell.x, p_cell.y);
	if (p_type == FLUID_CELL && s[idx] != 0.0) {
		_spawn_cell_particles(p_cell.x * cell_size, p_cell.y * cell_size);
		cell_types[idx] = FLUID_CELL;
	}
}

FLIPWaterSimulation::CellType FLIPWaterSimulation::get_cell(const Vector2i &p_cell) const {
	if (!_is_cell_in_bounds(p_cell.x, p_cell.y)) {
		return AIR_CELL;
	}
	return static_cast<CellType>(cell_types[p_cell.x + (f_num_x * p_cell.y)]);
}

void FLIPWaterSimulation::fill_cells(const Rect2i &p_rect, CellType p_type) {
	Rect2i rect = p_rect.intersection(Rect2i(0, 0, f_num_x, f_num_y));
	for (int j = rect.position.y; j < rect.get_end().y; j++) {
		for (int i = rect.position.x; i < rect.get_end().x; i++) {
			set_cell(Vector2i(i, j), p_type);
		}
	}
}

void FLIPWaterSimulation::clear() {
	p_positions.clear();
	p_velocities.clear();
	std::fill(user_solid.begin(), user_solid.end(), 0);
	std::fill(u.begin(), u.end(), 0.0);
	std::fill(v.begin(), v.end(), 0.0);
	std::fill(pressure.begin(), pressure.end(), 0.0);
	std::fill(cell_types.begin(), cell_types.end(), AIR_CELL);
	_update_solids();
	_update_render();
}

void FLIPWaterSimulation::shift(const Vector2i &p_cells) {
	if (p_cells == Vector2i()) {
		return;
	}
	int dx = p_cells.x;
	int dy = p_cells.y;
	int nx = f_num_x;
	int ny = f_num_y;

	shift_grid<uint8_t>(user_solid, nx, ny, dx, dy, 0);
	shift_grid<double>(u, nx, ny, dx, dy, 0.0);
	shift_grid<double>(v, nx, ny, dx, dy, 0.0);
	shift_grid<double>(prev_u, nx, ny, dx, dy, 0.0);
	shift_grid<double>(prev_v, nx, ny, dx, dy, 0.0);
	shift_grid<double>(pressure, nx, ny, dx, dy, 0.0);
	shift_grid<double>(p_densities, nx, ny, dx, dy, 0.0);
	shift_grid<int>(cell_types, nx, ny, dx, dy, AIR_CELL);
	// Collider layers are rebuilt from the physics world at their new position.
	std::fill(static_solid.begin(), static_solid.end(), 0);
	std::fill(dynamic_solid.begin(), dynamic_solid.end(), 0);
	std::fill(solid_velocity.begin(), solid_velocity.end(), Vector2());
	dynamic_cells.clear();
	static_colliders_dirty = true;

	Vector2 offset = Vector2(p_cells) * cell_size;
	for (size_t i = 0; i < p_positions.size(); i++) {
		p_positions[i] -= offset;
	}
	_remove_particles_outside_domain();
	set_position(get_position() + get_transform().basis_xform(offset));
	_update_solids();

	if (dx != 0) {
		int w = std::min(std::abs(dx), nx);
		emit_signal("region_exposed", Rect2i(dx > 0 ? nx - w : 0, 0, w, ny));
	}
	if (dy != 0) {
		int h = std::min(std::abs(dy), ny);
		emit_signal("region_exposed", Rect2i(0, dy > 0 ? ny - h : 0, nx, h));
	}
}

void FLIPWaterSimulation::_follow() {
	if (follow_target.is_empty()) {
		return;
	}
	Node2D *target = Object::cast_to<Node2D>(get_node_or_null(follow_target));
	if (!target) {
		return;
	}
	Vector2 diff = (to_local(target->get_global_position()) - get_domain_rect().get_center()) / cell_size;
	Vector2i cells(static_cast<int>(diff.x), static_cast<int>(diff.y));
	if (cells != Vector2i()) {
		shift(cells);
	}
}

// Particles

bool FLIPWaterSimulation::_add_particle_internal(const Vector2 &pos, const Vector2 &vel) {
	if (static_cast<int>(p_positions.size()) >= max_particles) {
		return false;
	}
	p_positions.push_back(pos);
	p_velocities.push_back(vel);
	return true;
}

void FLIPWaterSimulation::_remove_particles_outside_domain() {
	Rect2 domain = get_domain_rect();
	size_t kept = 0;
	for (size_t i = 0; i < p_positions.size(); i++) {
		if (domain.has_point(p_positions[i])) {
			p_positions[kept] = p_positions[i];
			p_velocities[kept] = p_velocities[i];
			kept++;
		}
	}
	p_positions.resize(kept);
	p_velocities.resize(kept);
}

int FLIPWaterSimulation::get_particle_count() const {
	return static_cast<int>(p_positions.size());
}

PackedVector2Array FLIPWaterSimulation::get_particle_positions() const {
	PackedVector2Array out;
	out.resize(p_positions.size());
	std::copy(p_positions.begin(), p_positions.end(), out.ptrw());
	return out;
}

PackedVector2Array FLIPWaterSimulation::get_particle_velocities() const {
	PackedVector2Array out;
	out.resize(p_velocities.size());
	std::copy(p_velocities.begin(), p_velocities.end(), out.ptrw());
	return out;
}

PackedInt32Array FLIPWaterSimulation::get_cell_types() const {
	PackedInt32Array out;
	out.resize(cell_types.size());
	std::copy(cell_types.begin(), cell_types.end(), out.ptrw());
	return out;
}

PackedFloat32Array FLIPWaterSimulation::get_pressure() const {
	PackedFloat32Array out;
	out.resize(pressure.size());
	std::copy(pressure.begin(), pressure.end(), out.ptrw());
	return out;
}

bool FLIPWaterSimulation::add_particle(const Vector2 &p_pos, const Vector2 &p_velocity) {
	return _add_particle_internal(p_pos, p_velocity);
}

void FLIPWaterSimulation::apply_impulse(const Vector2 &p_center, double p_radius, const Vector2 &p_velocity) {
	double r2 = p_radius * p_radius;
	for (size_t i = 0; i < p_positions.size(); i++) {
		if (p_positions[i].distance_squared_to(p_center) <= r2) {
			p_velocities[i] += p_velocity;
		}
	}
}

int FLIPWaterSimulation::remove_particles_in_circle(const Vector2 &p_center, double p_radius) {
	double r2 = p_radius * p_radius;
	size_t kept = 0;
	for (size_t i = 0; i < p_positions.size(); i++) {
		if (p_positions[i].distance_squared_to(p_center) > r2) {
			p_positions[kept] = p_positions[i];
			p_velocities[kept] = p_velocities[i];
			kept++;
		}
	}
	int removed = static_cast<int>(p_positions.size() - kept);
	p_positions.resize(kept);
	p_velocities.resize(kept);
	return removed;
}

int FLIPWaterSimulation::count_particles_in_circle(const Vector2 &p_center, double p_radius) const {
	double r2 = p_radius * p_radius;
	int count = 0;
	for (size_t i = 0; i < p_positions.size(); i++) {
		if (p_positions[i].distance_squared_to(p_center) <= r2) {
			count++;
		}
	}
	return count;
}

int FLIPWaterSimulation::count_particles_in_rect(const Rect2 &p_rect) const {
	Rect2 rect = p_rect.abs();
	int count = 0;
	for (size_t i = 0; i < p_positions.size(); i++) {
		if (rect.has_point(p_positions[i])) {
			count++;
		}
	}
	return count;
}

// Physics colliders

void FLIPWaterSimulation::refresh_colliders() {
	static_colliders_dirty = true;
}

void FLIPWaterSimulation::_refresh_static_colliders() {
	static_colliders_dirty = false;
	std::fill(static_solid.begin(), static_solid.end(), 0);
	static_bodies.clear();
	static_body_rids.clear();
	if (!is_inside_tree() || collision_mask == 0) {
		return;
	}
	PhysicsDirectSpaceState2D *space = get_world_2d()->get_direct_space_state();
	if (!space) {
		return;
	}
	PhysicsServer2D *ps = PhysicsServer2D::get_singleton();
	if (!query_shape.is_valid()) {
		query_shape = ps->rectangle_shape_create();
		query_params.instantiate();
		query_params->set_shape_rid(query_shape);
		query_params->set_collide_with_areas(false);
		query_params->set_collide_with_bodies(true);
	}
	Rect2 domain = get_domain_rect();
	Transform2D sim_xform = get_global_transform();
	ps->shape_set_data(query_shape, domain.size * 0.5);
	query_params->set_transform(sim_xform * Transform2D(0.0, domain.get_center()));
	query_params->set_collision_mask(collision_mask);

	// Results are capped, so keep querying while excluding everything already seen.
	TypedArray<RID> exclude;
	std::unordered_set<uint64_t> seen_dynamic;
	for (int pass = 0; pass < 64; pass++) {
		query_params->set_exclude(exclude);
		TypedArray<Dictionary> hits = space->intersect_shape(query_params, MAX_QUERY_RESULTS);
		bool capped = hits.size() >= MAX_QUERY_RESULTS;
		for (int k = 0; k < hits.size(); k++) {
			Dictionary hit = hits[k];
			RID rid = hit["rid"];
			int shape = hit["shape"];
			uint64_t id = rid.get_id();
			if (ps->body_get_mode(rid) != PhysicsServer2D::BODY_MODE_STATIC) {
				if (seen_dynamic.insert(id).second) {
					exclude.push_back(rid);
				}
				continue;
			}
			Transform2D body_xform = ps->body_get_state(rid, PhysicsServer2D::BODY_STATE_TRANSFORM);
			bool is_new = static_bodies.insert(id).second;
			if (is_new) {
				static_body_rids.push_back(rid);
				exclude.push_back(rid);
			}
			if (!capped) {
				_rasterize_shape(rid, shape, false, sim_xform, body_xform, 0.0);
			} else if (is_new) {
				// This body is excluded from the next pass, so take all of its shapes now.
				for (int n = 0; n < ps->body_get_shape_count(rid); n++) {
					_rasterize_shape(rid, n, false, sim_xform, body_xform, 0.0);
				}
			}
		}
		if (!capped) {
			break;
		}
	}
	query_params->set_exclude(static_body_rids);
}

void FLIPWaterSimulation::_update_dynamic_colliders(double delta) {
	for (int idx : dynamic_cells) {
		dynamic_solid[idx] = 0;
		solid_velocity[idx] = Vector2();
	}
	dynamic_cells.clear();
	if (!is_inside_tree() || collision_mask == 0 || !query_params.is_valid()) {
		return;
	}
	PhysicsDirectSpaceState2D *space = get_world_2d()->get_direct_space_state();
	if (!space) {
		return;
	}
	PhysicsServer2D *ps = PhysicsServer2D::get_singleton();
	Rect2 domain = get_domain_rect();
	Transform2D sim_xform = get_global_transform();
	ps->shape_set_data(query_shape, domain.size * 0.5);
	query_params->set_transform(sim_xform * Transform2D(0.0, domain.get_center()));
	query_params->set_collision_mask(collision_mask);

	TypedArray<Dictionary> hits = space->intersect_shape(query_params, MAX_QUERY_RESULTS);
	std::unordered_map<uint64_t, Transform2D> current_transforms;
	bool new_static = false;
	for (int k = 0; k < hits.size(); k++) {
		Dictionary hit = hits[k];
		RID rid = hit["rid"];
		int shape = hit["shape"];
		uint64_t id = rid.get_id();
		Transform2D body_xform = ps->body_get_state(rid, PhysicsServer2D::BODY_STATE_TRANSFORM);
		if (ps->body_get_mode(rid) == PhysicsServer2D::BODY_MODE_STATIC) {
			// A static body added since the last refresh.
			if (static_bodies.insert(id).second) {
				static_body_rids.push_back(rid);
				new_static = true;
			}
			_rasterize_shape(rid, shape, false, sim_xform, body_xform, 0.0);
			continue;
		}
		current_transforms[id] = body_xform;
		_rasterize_shape(rid, shape, true, sim_xform, body_xform, delta);
	}
	body_prev_transforms.swap(current_transforms);
	if (new_static) {
		query_params->set_exclude(static_body_rids);
	}
}

void FLIPWaterSimulation::_rasterize_shape(const RID &body, int shape_idx, bool dynamic, const Transform2D &sim_xform, const Transform2D &body_xform, double delta) {
	PhysicsServer2D *ps = PhysicsServer2D::get_singleton();
	RID shape = ps->body_get_shape(body, shape_idx);
	Transform2D sim_inv = sim_xform.affine_inverse();
	Transform2D local_xf = sim_inv * body_xform * ps->body_get_shape_transform(body, shape_idx);
	Transform2D to_shape = local_xf.affine_inverse();
	Variant data = ps->shape_get_data(shape);
	double h = cell_size;
	double shape_scale = 0.5 * (to_shape.basis_xform(Vector2(1, 0)).length() + to_shape.basis_xform(Vector2(0, 1)).length());

	// Cells count as solid when their centre is within `margin` of the shape. Thin shapes (segments)
	// get a wider margin so they always block at least one cell.
	double area_margin = 0.25 * h * shape_scale;
	double line_margin = 0.5 * h * shape_scale;
	double margin = area_margin;
	Rect2 bounds;
	bool unbounded = false;
	std::vector<Vector2> points;

	PhysicsServer2D::ShapeType type = ps->shape_get_type(shape);
	switch (type) {
		case PhysicsServer2D::SHAPE_RECTANGLE: {
			Vector2 he = data;
			bounds = Rect2(-he, he * 2.0);
		} break;
		case PhysicsServer2D::SHAPE_CIRCLE: {
			double r = data;
			bounds = Rect2(-r, -r, 2.0 * r, 2.0 * r);
		} break;
		case PhysicsServer2D::SHAPE_CAPSULE: {
			Vector2 hr = data; // (height, radius)
			bounds = Rect2(-hr.y, -hr.x * 0.5, 2.0 * hr.y, hr.x);
		} break;
		case PhysicsServer2D::SHAPE_CONVEX_POLYGON:
		case PhysicsServer2D::SHAPE_CONCAVE_POLYGON: {
			PackedVector2Array pts = data;
			if (pts.is_empty()) {
				return;
			}
			points.assign(pts.ptr(), pts.ptr() + pts.size());
			bounds = Rect2(points[0], Vector2());
			for (const Vector2 &p : points) {
				bounds.expand_to(p);
			}
			if (type == PhysicsServer2D::SHAPE_CONCAVE_POLYGON) {
				margin = line_margin;
			}
		} break;
		case PhysicsServer2D::SHAPE_SEGMENT: {
			Rect2 seg = data; // position = a, size = b
			points = { seg.position, seg.size };
			bounds = Rect2(seg.position, Vector2());
			bounds.expand_to(seg.size);
			margin = line_margin;
		} break;
		case PhysicsServer2D::SHAPE_WORLD_BOUNDARY:
			unbounded = true;
			break;
		default:
			return;
	}

	// Signed distance from a point in shape space to the shape.
	auto sdf = [&](const Vector2 &p) -> double {
		switch (type) {
			case PhysicsServer2D::SHAPE_RECTANGLE: {
				Vector2 he = data;
				Vector2 q = p.abs() - he;
				return Vector2(std::max(q.x, real_t(0)), std::max(q.y, real_t(0))).length() + std::min(std::max(q.x, q.y), real_t(0));
			}
			case PhysicsServer2D::SHAPE_CIRCLE:
				return p.length() - static_cast<double>(data);
			case PhysicsServer2D::SHAPE_CAPSULE: {
				Vector2 hr = data;
				double half = std::max(hr.x * 0.5 - hr.y, 0.0);
				return Vector2(p.x, p.y - std::clamp(static_cast<double>(p.y), -half, half)).length() - hr.y;
			}
			case PhysicsServer2D::SHAPE_CONVEX_POLYGON: {
				double dist = std::numeric_limits<double>::max();
				bool inside = false;
				size_t n = points.size();
				for (size_t a = 0, b = n - 1; a < n; b = a++) {
					const Vector2 &pa = points[a];
					const Vector2 &pb = points[b];
					if (((pa.y > p.y) != (pb.y > p.y)) && (p.x < (pb.x - pa.x) * (p.y - pa.y) / (pb.y - pa.y) + pa.x)) {
						inside = !inside;
					}
					dist = std::min(dist, segment_distance(p, pa, pb));
				}
				return inside ? -dist : dist;
			}
			case PhysicsServer2D::SHAPE_CONCAVE_POLYGON:
			case PhysicsServer2D::SHAPE_SEGMENT: {
				double dist = std::numeric_limits<double>::max();
				for (size_t a = 0; a + 1 < points.size(); a += 2) {
					dist = std::min(dist, segment_distance(p, points[a], points[a + 1]));
				}
				return dist;
			}
			case PhysicsServer2D::SHAPE_WORLD_BOUNDARY: {
				Array arr = data;
				Vector2 normal = arr[0];
				double d = arr[1];
				return normal.dot(p) - d;
			}
			default:
				return 1.0;
		}
	};

	int i0 = 0, j0 = 0, i1 = f_num_x - 1, j1 = f_num_y - 1;
	if (!unbounded) {
		Rect2 local_bounds = local_xf.xform(bounds).grow(0.5 * h);
		i0 = std::max(static_cast<int>(std::floor(local_bounds.position.x / h)), 0);
		j0 = std::max(static_cast<int>(std::floor(local_bounds.position.y / h)), 0);
		i1 = std::min(static_cast<int>(std::floor(local_bounds.get_end().x / h)), f_num_x - 1);
		j1 = std::min(static_cast<int>(std::floor(local_bounds.get_end().y / h)), f_num_y - 1);
	}

	const Transform2D *prev_xform = nullptr;
	Transform2D body_inv;
	if (dynamic && delta > 0.0) {
		auto it = body_prev_transforms.find(body.get_id());
		if (it != body_prev_transforms.end()) {
			prev_xform = &it->second;
			body_inv = body_xform.affine_inverse();
		}
	}

	for (int j = j0; j <= j1; j++) {
		for (int i = i0; i <= i1; i++) {
			Vector2 center((i + 0.5) * h, (j + 0.5) * h);
			if (sdf(to_shape.xform(center)) > margin) {
				continue;
			}
			int idx = i + (f_num_x * j);
			if (!dynamic) {
				static_solid[idx] = 1;
				continue;
			}
			if (!dynamic_solid[idx]) {
				dynamic_solid[idx] = 1;
				dynamic_cells.push_back(idx);
			}
			if (prev_xform) {
				// Velocity of the body's material point now at this cell, from its motion since the last frame.
				Vector2 global = sim_xform.xform(center);
				Vector2 global_prev = prev_xform->xform(body_inv.xform(global));
				solid_velocity[idx] = sim_inv.basis_xform((global - global_prev) / delta);
			}
		}
	}
}

// Effectors

void FLIPWaterSimulation::_gather_effectors() {
	frame_effectors.clear();
	frame_effector_ids.clear();
	if (!is_inside_tree()) {
		return;
	}
	TypedArray<Node> nodes = get_tree()->get_nodes_in_group(FLIPEffector::GROUP);
	for (int k = 0; k < nodes.size(); k++) {
		FLIPEffector *effector = Object::cast_to<FLIPEffector>(nodes[k]);
		if (effector && effector->is_enabled() && effector->get_simulation_node() == this) {
			frame_effectors.push_back(effector);
			frame_effector_ids.push_back(effector->get_instance_id());
		}
	}
}

void FLIPWaterSimulation::_apply_effectors(double dt) {
	if (frame_effectors.empty()) {
		return;
	}
	Transform2D sim_inv = get_global_transform().affine_inverse();
	for (FLIPEffector *effector : frame_effectors) {
		Transform2D xf = sim_inv * effector->get_global_transform();
		if (FLIPEmitter *emitter = Object::cast_to<FLIPEmitter>(effector)) {
			_apply_emitter(emitter, xf, dt);
		} else if (FLIPForce *force = Object::cast_to<FLIPForce>(effector)) {
			_apply_force(force, xf, dt);
		} else if (FLIPDrain *drain = Object::cast_to<FLIPDrain>(effector)) {
			drain->add_drained(_apply_drain(drain, xf));
		}
	}
}

void FLIPWaterSimulation::_apply_emitter(FLIPEmitter *emitter, const Transform2D &xf, double dt) {
	int count = emitter->take_spawn_count(dt);
	Vector2 dir = xf.basis_xform(Vector2(1, 0)).normalized();
	std::uniform_real_distribution<double> unit(-0.5, 0.5);
	for (int k = 0; k < count; k++) {
		// Jitter along the emission direction too, so particles never spawn on top of each other.
		Vector2 pos = xf.xform(Vector2(unit(rng) * p_radius, unit(rng) * emitter->get_width()));
		Vector2i cell = local_to_cell(pos);
		if (!_is_cell_in_bounds(cell.x, cell.y) || s[cell.x + (f_num_x * cell.y)] == 0.0) {
			continue;
		}
		Vector2 vel = dir.rotated(unit(rng) * emitter->get_spread()) * emitter->get_speed();
		if (!_add_particle_internal(pos, vel)) {
			break;
		}
	}
}

void FLIPWaterSimulation::_apply_force(FLIPForce *force, const Transform2D &xf, double dt) {
	Vector2 center = xf.get_origin();
	Vector2 dir = xf.basis_xform(Vector2(1, 0));
	double scale = dir.length();
	if (scale == 0.0) {
		return;
	}
	dir /= scale;
	double radius = force->get_radius() * scale;
	double r2 = radius * radius;
	bool full_circle = force->get_spread() >= Math_TAU - 0.001;
	double cos_half_spread = std::cos(force->get_spread() * 0.5);
	bool radial = force->get_mode() == FLIPForce::MODE_RADIAL;
	double strength = force->get_strength() * dt;

	for (size_t i = 0; i < p_positions.size(); i++) {
		Vector2 offset = p_positions[i] - center;
		double d2 = offset.length_squared();
		if (d2 > r2) {
			continue;
		}
		double d = std::sqrt(d2);
		if (!full_circle && (d == 0.0 || offset.dot(dir) / d < cos_half_spread)) {
			continue;
		}
		Vector2 accel_dir = radial ? (d > 0.0 ? offset / d : Vector2()) : dir;
		double factor = force->has_falloff() ? 1.0 - (d / radius) : 1.0;
		p_velocities[i] += accel_dir * strength * factor;
	}
}

int FLIPWaterSimulation::_apply_drain(FLIPDrain *drain, const Transform2D &xf) {
	double scale = xf.basis_xform(Vector2(1, 0)).length();
	return remove_particles_in_circle(xf.get_origin(), drain->get_radius() * scale);
}

void FLIPWaterSimulation::_update_sensors() {
	Transform2D sim_xform = get_global_transform();
	for (uint64_t id : frame_effector_ids) {
		// Re-resolve every time: a signal handler may have freed any of these nodes.
		Object *obj = ObjectDB::get_instance(id);
		if (FLIPSensor *sensor = Object::cast_to<FLIPSensor>(obj)) {
			if (!sensor->is_inside_tree()) {
				continue;
			}
			Transform2D to_sensor = sensor->get_global_transform().affine_inverse() * sim_xform;
			Vector2 half = sensor->get_size() * 0.5;
			int count = 0;
			for (size_t i = 0; i < p_positions.size(); i++) {
				Vector2 p = to_sensor.xform(p_positions[i]);
				if (std::abs(p.x) <= half.x && std::abs(p.y) <= half.y) {
					count++;
				}
			}
			sensor->set_particle_count(count);
		} else if (FLIPDrain *drain = Object::cast_to<FLIPDrain>(obj)) {
			drain->flush_drained();
		}
	}
}

// Simulation

void FLIPWaterSimulation::update(double delta) {
	if (delta <= 0.0) {
		return;
	}

	if (is_inside_tree() && collision_mask != 0) {
		if (static_colliders_dirty) {
			_refresh_static_colliders();
		}
		_update_dynamic_colliders(delta);
	}
	_update_solids();
	_gather_effectors();

	double dt = delta / substeps;
	for (int i = 0; i < substeps; i++) {
		_step(dt);
	}

	frame_effectors.clear();
	_update_sensors();
	frame_effector_ids.clear();
	_update_render();
	if (debug_draw) {
		queue_redraw();
	}
}

void FLIPWaterSimulation::_step(double dt) {
	_apply_effectors(dt);
	_integrate_particles(dt);
	_push_particles_apart();
	_handle_collisions();
	_particles_to_grid();
	_update_particle_density();
	_solve_incompressibility(dt);
	_run_parallel(static_cast<int>(p_positions.size()), 2048, &FLIPWaterSimulation::_grid_to_particles_range);
}

void FLIPWaterSimulation::_run_parallel(int count, int min_per_chunk, void (FLIPWaterSimulation::*fn)(int, int)) {
	if (count <= 0) {
		return;
	}
	int chunks = 1;
	if (use_threads) {
		chunks = std::min(OS::get_singleton()->get_processor_count(), count / std::max(min_per_chunk, 1));
	}
	if (chunks <= 1) {
		(this->*fn)(0, count);
		return;
	}
	job_fn = fn;
	job_count = count;
	job_chunks = chunks;
	WorkerThreadPool *pool = WorkerThreadPool::get_singleton();
	int64_t group = pool->add_group_task(callable_mp(this, &FLIPWaterSimulation::_run_job_chunk), chunks, chunks, true);
	pool->wait_for_group_task_completion(group);
}

void FLIPWaterSimulation::_run_job_chunk(uint32_t chunk) {
	int begin = static_cast<int>((static_cast<int64_t>(job_count) * chunk) / job_chunks);
	int end = static_cast<int>((static_cast<int64_t>(job_count) * (chunk + 1)) / job_chunks);
	(this->*job_fn)(begin, end);
}

void FLIPWaterSimulation::_integrate_particles(double dt) {
	// Never move more than one cell per step, so particles can't tunnel through walls.
	double limit = cell_size / dt;
	if (max_speed > 0.0) {
		limit = std::min(limit, max_speed);
	}

	for (size_t i = 0; i < p_positions.size(); i++) {
		p_velocities[i] += gravity * dt;
		double speed = p_velocities[i].length();
		if (speed > limit) {
			p_velocities[i] *= limit / speed;
		}
		p_positions[i] += p_velocities[i] * dt;
	}
}

void FLIPWaterSimulation::_handle_collisions() {
	double h = cell_size;
	double min_x = h + p_radius;
	double max_x = ((f_num_x - 1) * h) - p_radius;
	double min_y = h + p_radius;
	double max_y = ((f_num_y - 1) * h) - p_radius;

	// Compact surviving particles in place.
	size_t kept = 0;
	for (size_t i = 0; i < p_positions.size(); i++) {
		double x = p_positions[i].x;
		double y = p_positions[i].y;
		Vector2 vel = p_velocities[i];

		// Push particles out of solid cells (moving colliders, newly placed walls) into the nearest free cell.
		int ci = std::clamp(static_cast<int>(std::floor(x / h)), 0, f_num_x - 1);
		int cj = std::clamp(static_cast<int>(std::floor(y / h)), 0, f_num_y - 1);
		int cell = ci + (f_num_x * cj);
		if (s[cell] == 0.0) {
			double best_d2 = std::numeric_limits<double>::max();
			Vector2 best;
			for (int ring = 1; ring <= 2 && best_d2 == std::numeric_limits<double>::max(); ring++) {
				for (int nj = cj - ring; nj <= cj + ring; nj++) {
					for (int ni = ci - ring; ni <= ci + ring; ni++) {
						if (!_is_cell_in_bounds(ni, nj) || s[ni + (f_num_x * nj)] == 0.0) {
							continue;
						}
						Vector2 target(
								std::clamp(x, ni * h + p_radius, (ni + 1) * h - p_radius),
								std::clamp(y, nj * h + p_radius, (nj + 1) * h - p_radius));
						double d2 = target.distance_squared_to(Vector2(x, y));
						if (d2 < best_d2) {
							best_d2 = d2;
							best = target;
						}
					}
				}
			}
			if (best_d2 != std::numeric_limits<double>::max()) {
				x = best.x;
				y = best.y;
				vel = solid_velocity[cell];
			}
		}

		// Clamp to interior cells so the neighbour lookups stay in bounds for particles in the border cells.
		int pi = std::clamp(static_cast<int>(std::floor(x / h)), 1, f_num_x - 2);
		int pj = std::clamp(static_cast<int>(std::floor(y / h)), 1, f_num_y - 2);
		int index = pi + (f_num_x * pj);

		if (s[index - 1] == 0.0 && x - p_radius < pi * h) {
			x = (pi * h) + p_radius;
			vel.x = solid_velocity[index - 1].x;
		} else if (s[index + 1] == 0.0 && x + p_radius > (pi + 1) * h) {
			x = ((pi + 1) * h) - p_radius;
			vel.x = solid_velocity[index + 1].x;
		}

		if (s[index - f_num_x] == 0.0 && y - p_radius < pj * h) {
			y = (pj * h) + p_radius;
			vel.y = solid_velocity[index - f_num_x].y;
		} else if (s[index + f_num_x] == 0.0 && y + p_radius > (pj + 1) * h) {
			y = ((pj + 1) * h) - p_radius;
			vel.y = solid_velocity[index + f_num_x].y;
		}

		// Walls are resolved first so only particles escaping through an open edge are removed.
		if (x < min_x || x > max_x || y < min_y || y > max_y) {
			continue;
		}

		p_positions[kept] = Vector2(x, y);
		p_velocities[kept] = vel;
		kept++;
	}

	p_positions.resize(kept);
	p_velocities.resize(kept);
}

void FLIPWaterSimulation::_push_particles_apart() {
	int num_particles = static_cast<int>(p_positions.size());
	if (num_particles == 0) {
		return;
	}

	std::fill(num_cell_particles.begin(), num_cell_particles.end(), 0);
	cell_particle_ids.resize(num_particles);
	p_push.resize(num_particles);
	int p_num_cells = p_num_x * p_num_y;

	for (int i = 0; i < num_particles; i++) {
		int xi = std::clamp(static_cast<int>(std::floor(p_positions[i].x * p_inv_spacing)), 0, p_num_x - 1);
		int yi = std::clamp(static_cast<int>(std::floor(p_positions[i].y * p_inv_spacing)), 0, p_num_y - 1);
		num_cell_particles[xi + (p_num_x * yi)]++;
	}

	int _first = 0;
	for (int i = 0; i < p_num_cells; i++) {
		_first += num_cell_particles[i];
		first_cell_particle[i] = _first;
	}
	first_cell_particle[p_num_cells] = _first;

	for (int i = 0; i < num_particles; i++) {
		int xi = std::clamp(static_cast<int>(std::floor(p_positions[i].x * p_inv_spacing)), 0, p_num_x - 1);
		int yi = std::clamp(static_cast<int>(std::floor(p_positions[i].y * p_inv_spacing)), 0, p_num_y - 1);
		int cell_nr = xi + (p_num_x * yi);
		first_cell_particle[cell_nr]--;
		cell_particle_ids[first_cell_particle[cell_nr]] = i;
	}

	// Jacobi-style: each particle computes its own correction from a read-only snapshot, so this parallelises.
	for (int iter = 0; iter < 2; iter++) {
		_run_parallel(num_particles, 1024, &FLIPWaterSimulation::_push_particles_range);
		for (int i = 0; i < num_particles; i++) {
			p_positions[i] += p_push[i];
		}
	}
}

void FLIPWaterSimulation::_push_particles_range(int begin, int end) {
	double min_dist = 2.0 * p_radius;
	double min_dist2 = min_dist * min_dist;

	for (int i = begin; i < end; i++) {
		double x = p_positions[i].x;
		double y = p_positions[i].y;
		Vector2 push;

		int pxi = static_cast<int>(std::floor(x * p_inv_spacing));
		int pyi = static_cast<int>(std::floor(y * p_inv_spacing));
		int x0 = std::max(pxi - 1, 0);
		int y0 = std::max(pyi - 1, 0);
		int x1 = std::min(pxi + 1, p_num_x - 1);
		int y1 = std::min(pyi + 1, p_num_y - 1);

		for (int xi = x0; xi <= x1; xi++) {
			for (int yi = y0; yi <= y1; yi++) {
				int cell_nr = xi + (p_num_x * yi);
				int first = first_cell_particle[cell_nr];
				int last = first_cell_particle[cell_nr + 1];

				for (int j = first; j < last; j++) {
					int id = cell_particle_ids[j];
					if (id == i) {
						continue;
					}

					double dx = p_positions[id].x - x;
					double dy = p_positions[id].y - y;
					double d2 = dx * dx + dy * dy;
					if (d2 > min_dist2 || d2 == 0.0) {
						continue;
					}

					double d = std::sqrt(d2);
					double _s = 0.5 * (min_dist - d) / d;
					push.x -= dx * _s;
					push.y -= dy * _s;
				}
			}
		}
		p_push[i] = push;
	}
}

void FLIPWaterSimulation::_particles_to_grid() {
	double h = cell_size;
	double h1 = f_inv_spacing;
	double h2 = 0.5 * h;

	prev_u = u;
	prev_v = v;
	std::fill(du.begin(), du.end(), 0.0);
	std::fill(dv.begin(), dv.end(), 0.0);
	std::fill(u.begin(), u.end(), 0.0);
	std::fill(v.begin(), v.end(), 0.0);

	for (int i = 0; i < f_num_cells; i++) {
		cell_types[i] = s[i] == 0.0 ? SOLID_CELL : AIR_CELL;
	}

	for (size_t i = 0; i < p_positions.size(); i++) {
		int xi = std::clamp(static_cast<int>(std::floor(p_positions[i].x * h1)), 0, f_num_x - 1);
		int yi = std::clamp(static_cast<int>(std::floor(p_positions[i].y * h1)), 0, f_num_y - 1);
		int cell_nr = xi + (f_num_x * yi);
		if (cell_types[cell_nr] == AIR_CELL) {
			cell_types[cell_nr] = FLUID_CELL;
		}
	}

	for (int component = 0; component < 2; component++) {
		double dx = component == 0 ? 0.0 : h2;
		double dy = component == 0 ? h2 : 0.0;
		std::vector<double> &f = component == 0 ? u : v;
		std::vector<double> &d = component == 0 ? du : dv;

		for (size_t i = 0; i < p_positions.size(); i++) {
			double x = std::clamp(static_cast<double>(p_positions[i].x), h, f_num_x * h);
			double y = std::clamp(static_cast<double>(p_positions[i].y), h, f_num_y * h);

			int x0 = std::min(static_cast<int>(std::floor((x - dx) * h1)), f_num_x - 1);
			double tx = ((x - dx) - x0 * h) * h1;
			int x1 = std::min(x0 + 1, f_num_x - 1);

			int y0 = std::min(static_cast<int>(std::floor((y - dy) * h1)), f_num_y - 1);
			double ty = ((y - dy) - y0 * h) * h1;
			int y1 = std::min(y0 + 1, f_num_y - 1);

			double sx = 1.0 - tx;
			double sy = 1.0 - ty;
			double w[4] = { sx * sy, tx * sy, tx * ty, sx * ty };
			int nr[4] = { x0 + (f_num_x * y0), x1 + (f_num_x * y0), x1 + (f_num_x * y1), x0 + (f_num_x * y1) };
			double vel = component == 0 ? p_velocities[i].x : p_velocities[i].y;

			for (int k = 0; k < 4; k++) {
				if (cell_types[nr[k]] != SOLID_CELL) {
					f[nr[k]] += vel * w[k];
					d[nr[k]] += w[k];
				}
			}
		}

		for (size_t i = 0; i < f.size(); i++) {
			if (d[i] > 0.0) {
				f[i] /= d[i];
			}
		}
	}

	// Faces touching a solid take the solid's velocity (zero for static walls).
	for (int j = 0; j < f_num_y; j++) {
		for (int i = 0; i < f_num_x; i++) {
			int idx = i + (f_num_x * j);
			bool solid = cell_types[idx] == SOLID_CELL;
			bool left_solid = i > 0 && cell_types[idx - 1] == SOLID_CELL;
			bool top_solid = j > 0 && cell_types[idx - f_num_x] == SOLID_CELL;
			if (solid || left_solid) {
				u[idx] = solid && left_solid ? 0.5 * (solid_velocity[idx].x + solid_velocity[idx - 1].x) : (solid ? solid_velocity[idx].x : solid_velocity[idx - 1].x);
			}
			if (solid || top_solid) {
				v[idx] = solid && top_solid ? 0.5 * (solid_velocity[idx].y + solid_velocity[idx - f_num_x].y) : (solid ? solid_velocity[idx].y : solid_velocity[idx - f_num_x].y);
			}
		}
	}
}

void FLIPWaterSimulation::_grid_to_particles_range(int begin, int end) {
	double h = cell_size;
	double h1 = f_inv_spacing;
	double h2 = 0.5 * h;

	for (int i = begin; i < end; i++) {
		for (int component = 0; component < 2; component++) {
			double dx = component == 0 ? 0.0 : h2;
			double dy = component == 0 ? h2 : 0.0;
			const std::vector<double> &f = component == 0 ? u : v;
			const std::vector<double> &prev_f = component == 0 ? prev_u : prev_v;

			double x = std::clamp(static_cast<double>(p_positions[i].x), h, f_num_x * h);
			double y = std::clamp(static_cast<double>(p_positions[i].y), h, f_num_y * h);

			int x0 = std::min(static_cast<int>(std::floor((x - dx) * h1)), f_num_x - 1);
			double tx = ((x - dx) - x0 * h) * h1;
			int x1 = std::min(x0 + 1, f_num_x - 1);

			int y0 = std::min(static_cast<int>(std::floor((y - dy) * h1)), f_num_y - 1);
			double ty = ((y - dy) - y0 * h) * h1;
			int y1 = std::min(y0 + 1, f_num_y - 1);

			double sx = 1.0 - tx;
			double sy = 1.0 - ty;
			double w[4] = { sx * sy, tx * sy, tx * ty, sx * ty };
			int nr[4] = { x0 + (f_num_x * y0), x1 + (f_num_x * y0), x1 + (f_num_x * y1), x0 + (f_num_x * y1) };
			int offset = component == 1 ? f_num_x : 1;

			double total = 0.0;
			double pic_v = 0.0;
			double corr = 0.0;
			for (int k = 0; k < 4; k++) {
				if (cell_types[nr[k]] != SOLID_CELL && cell_types[nr[k] - offset] != SOLID_CELL) {
					total += w[k];
					pic_v += w[k] * f[nr[k]];
					corr += w[k] * (f[nr[k]] - prev_f[nr[k]]);
				}
			}

			if (total > 0.0) {
				double vel = component == 0 ? p_velocities[i].x : p_velocities[i].y;
				double new_v = ((1.0 - flip_ratio) * (pic_v / total)) + flip_ratio * (vel + corr / total);
				if (component == 0) {
					p_velocities[i].x = new_v;
				} else {
					p_velocities[i].y = new_v;
				}
			}
		}
	}
}

void FLIPWaterSimulation::_update_particle_density() {
	std::fill(p_densities.begin(), p_densities.end(), 0.0);

	for (size_t i = 0; i < p_positions.size(); i++) {
		double x = std::clamp(static_cast<double>(p_positions[i].x), cell_size, (f_num_x - 1) * cell_size);
		double y = std::clamp(static_cast<double>(p_positions[i].y), cell_size, (f_num_y - 1) * cell_size);

		int x0 = static_cast<int>(std::floor((x - (0.5 * cell_size)) * f_inv_spacing));
		double tx = ((x - (0.5 * cell_size)) - x0 * cell_size) * f_inv_spacing;
		int x1 = std::min(x0 + 1, f_num_x - 1);

		int y0 = static_cast<int>(std::floor((y - (0.5 * cell_size)) * f_inv_spacing));
		double ty = ((y - (0.5 * cell_size)) - y0 * cell_size) * f_inv_spacing;
		int y1 = std::min(y0 + 1, f_num_y - 1);

		double sx = 1.0 - tx;
		double sy = 1.0 - ty;

		p_densities[x0 + (f_num_x * y0)] += sx * sy;
		p_densities[x1 + (f_num_x * y0)] += tx * sy;
		p_densities[x1 + (f_num_x * y1)] += tx * ty;
		p_densities[x0 + (f_num_x * y1)] += sx * ty;
	}

	if (p_rest_density == 0.0) {
		double sum = 0.0;
		double num_fluid_cells = 0.0;

		for (int i = 0; i < f_num_cells; i++) {
			if (cell_types[i] == FLUID_CELL) {
				sum += p_densities[i];
				num_fluid_cells++;
			}
		}

		if (num_fluid_cells > 0.0) {
			p_rest_density = sum / num_fluid_cells;
		}
	}
}

void FLIPWaterSimulation::_solve_incompressibility(double dt) {
	std::fill(pressure.begin(), pressure.end(), 0.0);
	prev_u = u;
	prev_v = v;
	solve_cp = DENSITY * cell_size / dt;

	// Red-black Gauss-Seidel: cells of one colour share no faces, so each half-sweep can run in parallel.
	// Only worth threading on large grids, so each chunk must cover at least ~8k cells.
	int min_rows = std::max(1, 8192 / f_num_x);
	for (int iter = 0; iter < solver_iterations; iter++) {
		for (int color = 0; color < 2; color++) {
			solve_color = color;
			_run_parallel(f_num_y - 2, min_rows, &FLIPWaterSimulation::_solve_rows_range);
		}
	}
}

void FLIPWaterSimulation::_solve_rows_range(int begin, int end) {
	double drift = DRIFT_COMPENSATION * cell_size;

	for (int j = begin + 1; j < end + 1; j++) {
		int i0 = ((1 + j) % 2 == solve_color) ? 1 : 2;
		for (int i = i0; i < f_num_x - 1; i += 2) {
			int center = i + (f_num_x * j);
			if (cell_types[center] != FLUID_CELL) {
				continue;
			}

			int left = center - 1;
			int right = center + 1;
			int top = center - f_num_x;
			int bottom = center + f_num_x;

			double sx0 = s[left];
			double sx1 = s[right];
			double sy0 = s[top];
			double sy1 = s[bottom];
			double _s = sx0 + sx1 + sy0 + sy1;
			if (_s == 0.0) {
				continue;
			}

			double div = u[right] - u[center] + v[bottom] - v[center];
			if (p_rest_density > 0.0) {
				double compression = p_densities[center] - p_rest_density;
				if (compression > 0.0) {
					div -= drift * compression;
				}
			}

			double _p = -div / _s;
			_p *= OVER_RELAXATION;
			pressure[center] += solve_cp * _p;

			u[center] -= sx0 * _p;
			u[right] += sx1 * _p;
			v[center] -= sy0 * _p;
			v[bottom] += sy1 * _p;
		}
	}
}

// Rendering

void FLIPWaterSimulation::_ensure_render_resources() {
	if (render_item.is_valid()) {
		return;
	}
	RenderingServer *rs = RenderingServer::get_singleton();
	render_item = rs->canvas_item_create();
	render_multimesh = rs->multimesh_create();
	render_quad.instantiate();
	rs->multimesh_set_mesh(render_multimesh, render_quad->get_rid());
	rs->canvas_item_add_multimesh(render_item, render_multimesh);

	particle_material = make_material(PARTICLE_SHADER);
	splat_material = make_material(SPLAT_SHADER);
	surface_material = make_material(SURFACE_SHADER);

	surface_viewport = memnew(SubViewport);
	surface_viewport->set_transparent_background(true);
	surface_viewport->set_use_hdr_2d(true);
	surface_viewport->set_disable_3d(true);
	surface_viewport->set_update_mode(SubViewport::UPDATE_DISABLED);
	add_child(surface_viewport, false, INTERNAL_MODE_BACK);

	surface_root = memnew(Node2D);
	surface_viewport->add_child(surface_root);

	surface_sprite = memnew(Sprite2D);
	surface_sprite->set_centered(false);
	surface_sprite->set_texture(surface_viewport->get_texture());
	surface_sprite->set_material(surface_material);
	surface_sprite->set_texture_filter(TEXTURE_FILTER_LINEAR);
	surface_sprite->set_visible(false);
	add_child(surface_sprite, false, INTERNAL_MODE_BACK);
}

void FLIPWaterSimulation::_free_render_resources() {
	if (!render_item.is_valid()) {
		return;
	}
	RenderingServer *rs = RenderingServer::get_singleton();
	rs->free_rid(render_item);
	rs->free_rid(render_multimesh);
	render_item = RID();
	render_multimesh = RID();
}

void FLIPWaterSimulation::_update_render() {
	if (!render_item.is_valid()) {
		return;
	}
	RenderingServer *rs = RenderingServer::get_singleton();
	bool surface = render_mode == RENDER_SURFACE;
	bool visible = render_mode != RENDER_NONE;
	rs->canvas_item_set_visible(render_item, visible);
	surface_sprite->set_visible(surface);
	surface_viewport->set_update_mode(surface ? SubViewport::UPDATE_ALWAYS : SubViewport::UPDATE_DISABLED);
	if (!visible) {
		return;
	}

	rs->canvas_item_set_parent(render_item, surface ? surface_root->get_canvas_item() : get_canvas_item());
	rs->canvas_item_set_material(render_item, (surface ? splat_material : particle_material)->get_rid());
	if (surface) {
		Vector2 domain = get_domain_rect().size;
		Vector2i size(std::max(1, static_cast<int>(std::ceil(domain.x * surface_resolution))), std::max(1, static_cast<int>(std::ceil(domain.y * surface_resolution))));
		if (surface_viewport->get_size() != size) {
			surface_viewport->set_size(size);
		}
		surface_root->set_scale(Vector2(surface_resolution, surface_resolution));
		surface_sprite->set_scale(Vector2(1.0 / surface_resolution, 1.0 / surface_resolution));
		surface_material->set_shader_parameter("water_color", water_color);
		surface_material->set_shader_parameter("foam_color", foam_color);
	}

	int count = static_cast<int>(p_positions.size());
	if (count > render_capacity) {
		render_capacity = std::max({ count, render_capacity * 2, 256 });
		rs->multimesh_allocate_data(render_multimesh, render_capacity, RenderingServer::MULTIMESH_TRANSFORM_2D, true);
		render_buffer.resize(render_capacity * RENDER_FLOATS_PER_INSTANCE);
		render_buffer.fill(0.0f);
	}

	// Surface splats are three times the particle size so neighbours blend into one blob.
	float size = static_cast<float>((surface ? 6.0 : 2.0) * p_radius * particle_scale);
	float *w = render_buffer.ptrw();
	for (int i = 0; i < count; i++) {
		float foam = foam_speed > 0.0 ? static_cast<float>(std::min(p_velocities[i].length() / foam_speed, 1.0)) : 0.0f;
		Color c = surface ? Color(foam, 0, 0, 1) : water_color.lerp(foam_color, foam);
		float *o = w + i * RENDER_FLOATS_PER_INSTANCE;
		o[0] = size;
		o[1] = 0.0f;
		o[2] = 0.0f;
		o[3] = p_positions[i].x;
		o[4] = 0.0f;
		o[5] = size;
		o[6] = 0.0f;
		o[7] = p_positions[i].y;
		o[8] = c.r;
		o[9] = c.g;
		o[10] = c.b;
		o[11] = c.a;
	}
	rs->multimesh_set_buffer(render_multimesh, render_buffer);
	rs->multimesh_set_visible_instances(render_multimesh, count);
}
