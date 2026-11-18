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
	ClassDB::bind_method(D_METHOD("set_drift_compensation", "drift_compensation"), &FLIPWaterSimulation::set_drift_compensation);
	ClassDB::bind_method(D_METHOD("get_drift_compensation"), &FLIPWaterSimulation::get_drift_compensation);
	ClassDB::bind_method(D_METHOD("set_max_speed", "max_speed"), &FLIPWaterSimulation::set_max_speed);
	ClassDB::bind_method(D_METHOD("get_max_speed"), &FLIPWaterSimulation::get_max_speed);
	ClassDB::bind_method(D_METHOD("set_solver_iterations", "iterations"), &FLIPWaterSimulation::set_solver_iterations);
	ClassDB::bind_method(D_METHOD("get_solver_iterations"), &FLIPWaterSimulation::get_solver_iterations);
	ClassDB::bind_method(D_METHOD("set_particles_per_side", "count"), &FLIPWaterSimulation::set_particles_per_side);
	ClassDB::bind_method(D_METHOD("get_particles_per_side"), &FLIPWaterSimulation::get_particles_per_side);
	ClassDB::bind_method(D_METHOD("set_cohesion", "amount"), &FLIPWaterSimulation::set_cohesion);
	ClassDB::bind_method(D_METHOD("set_thin_runoff_ground", "enabled"), &FLIPWaterSimulation::set_thin_runoff_ground);
	ClassDB::bind_method(D_METHOD("get_thin_runoff_ground"), &FLIPWaterSimulation::get_thin_runoff_ground);
	ClassDB::bind_method(D_METHOD("get_cohesion"), &FLIPWaterSimulation::get_cohesion);
	ClassDB::bind_method(D_METHOD("set_viscosity", "amount"), &FLIPWaterSimulation::set_viscosity);
	ClassDB::bind_method(D_METHOD("get_viscosity"), &FLIPWaterSimulation::get_viscosity);
	ClassDB::bind_method(D_METHOD("set_separation_iterations", "iterations"), &FLIPWaterSimulation::set_separation_iterations);
	ClassDB::bind_method(D_METHOD("get_separation_iterations"), &FLIPWaterSimulation::get_separation_iterations);
	ClassDB::bind_method(D_METHOD("set_drift_pull", "enabled"), &FLIPWaterSimulation::set_drift_pull);
	ClassDB::bind_method(D_METHOD("is_drift_pull"), &FLIPWaterSimulation::is_drift_pull);
	ClassDB::bind_method(D_METHOD("set_spray_density", "fraction"), &FLIPWaterSimulation::set_spray_density);
	ClassDB::bind_method(D_METHOD("get_spray_density"), &FLIPWaterSimulation::get_spray_density);
	ClassDB::bind_method(D_METHOD("set_churn_low", "speed"), &FLIPWaterSimulation::set_churn_low);
	ClassDB::bind_method(D_METHOD("get_churn_low"), &FLIPWaterSimulation::get_churn_low);
	ClassDB::bind_method(D_METHOD("set_churn_high", "speed"), &FLIPWaterSimulation::set_churn_high);
	ClassDB::bind_method(D_METHOD("get_churn_high"), &FLIPWaterSimulation::get_churn_high);
	ClassDB::bind_method(D_METHOD("set_churn_fade", "rate"), &FLIPWaterSimulation::set_churn_fade);
	ClassDB::bind_method(D_METHOD("get_churn_fade"), &FLIPWaterSimulation::get_churn_fade);
	ClassDB::bind_method(D_METHOD("set_wall_drag", "rate"), &FLIPWaterSimulation::set_wall_drag);
	ClassDB::bind_method(D_METHOD("get_wall_drag"), &FLIPWaterSimulation::get_wall_drag);
	ClassDB::bind_method(D_METHOD("set_void_fill_cells", "cells"), &FLIPWaterSimulation::set_void_fill_cells);
	ClassDB::bind_method(D_METHOD("get_void_fill_cells"), &FLIPWaterSimulation::get_void_fill_cells);
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
	ClassDB::bind_method(D_METHOD("get_average_velocity_in_circle", "center", "radius"), &FLIPWaterSimulation::get_average_velocity_in_circle);
	ClassDB::bind_method(D_METHOD("apply_radial_impulse", "center", "radius", "strength"), &FLIPWaterSimulation::apply_radial_impulse);
	ClassDB::bind_method(D_METHOD("add_volume_source", "center", "radius", "rate"), &FLIPWaterSimulation::add_volume_source);
	ClassDB::bind_method(D_METHOD("build_surface", "origin", "samples", "spacing", "kernel", "radius", "band", "ghost_reach", "smooth", "level", "level_tol", "settle", "depth", "pocket"), &FLIPWaterSimulation::build_surface, DEFVAL(96.0), DEFVAL(900.0));
	ClassDB::bind_method(D_METHOD("sample_churned", "rect", "probes", "min_churn"), &FLIPWaterSimulation::sample_churned);
	ClassDB::bind_method(D_METHOD("sample_velocity", "position"), &FLIPWaterSimulation::sample_velocity);
	ClassDB::bind_method(D_METHOD("update_multimesh", "multimesh", "scale", "min_speed", "stretch", "max_stretch", "every", "min_churn", "mirror", "ground"), &FLIPWaterSimulation::update_multimesh, DEFVAL(1.0), DEFVAL(0.0), DEFVAL(0.0), DEFVAL(1.0), DEFVAL(1), DEFVAL(-1.0), DEFVAL(0.0), DEFVAL(0.0));

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
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "drift_compensation", PROPERTY_HINT_RANGE, "0,10,0.01,or_greater"), "set_drift_compensation", "get_drift_compensation");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "max_speed", PROPERTY_HINT_RANGE, "0,5000,1,or_greater,suffix:px/s"), "set_max_speed", "get_max_speed");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "solver_iterations", PROPERTY_HINT_RANGE, "1,200,1"), "set_solver_iterations", "get_solver_iterations");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "particles_per_side", PROPERTY_HINT_RANGE, "1,4,1"), "set_particles_per_side", "get_particles_per_side");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "separation_iterations", PROPERTY_HINT_RANGE, "0,16,1"), "set_separation_iterations", "get_separation_iterations");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "cohesion", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_cohesion", "get_cohesion");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "thin_runoff_ground"), "set_thin_runoff_ground", "get_thin_runoff_ground");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "viscosity", PROPERTY_HINT_RANGE, "0,1,0.01"), "set_viscosity", "get_viscosity");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "drift_pull"), "set_drift_pull", "is_drift_pull");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "wall_drag", PROPERTY_HINT_RANGE, "0,60,0.1,suffix:/s"), "set_wall_drag", "get_wall_drag");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "void_fill_cells", PROPERTY_HINT_RANGE, "0,1024,1"), "set_void_fill_cells", "get_void_fill_cells");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "spray_density", PROPERTY_HINT_RANGE, "0,0.9,0.01"), "set_spray_density", "get_spray_density");
	ADD_GROUP("Churn", "churn_");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "churn_low", PROPERTY_HINT_RANGE, "0,2000,1,suffix:px/s"), "set_churn_low", "get_churn_low");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "churn_high", PROPERTY_HINT_RANGE, "1,4000,1,suffix:px/s"), "set_churn_high", "get_churn_high");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "churn_fade", PROPERTY_HINT_RANGE, "0,20,0.01,suffix:/s"), "set_churn_fade", "get_churn_fade");
	ADD_GROUP("", "");
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
	drift_compensation = 5.0;
	max_speed = 0.0;
	solver_iterations = 30;
	particles_per_side = 2;
	separation_iterations = 2;
	drift_pull = true;
	cohesion = 0.0;
	thin_runoff_ground = true;
	viscosity = 0.0;
	spray_density = 0.0;
	churn_low = 120.0;
	churn_high = 400.0;
	churn_fade = 1.5;
	wall_drag = 0.0;
	void_fill_cells = 0;
	step_dt = 1.0 / 60.0;
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
	// Water is spawned as an n x n lattice per cell, and particles keep apart by exactly that spacing,
	// so freshly filled water is already at rest density. A radius that disagreed with the spawn
	// spacing let gravity squeeze water by about a fifth, and drift compensation (which only ever
	// pushes apart) then inflated it into a sponge.
	p_radius = 0.5 * cell_size / particles_per_side;
	// Each particle spreads a total weight of 1 over the cell centres around it, so evenly spread
	// water reads n * n everywhere. Measuring it from the first frame's water instead counted the
	// half-empty surface cells and read low.
	p_rest_density = particles_per_side * particles_per_side;

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
	drift_source.assign(f_num_cells, 0.0);
	thin_push.assign(f_num_cells, 0.0f);
	region_ids.assign(f_num_cells, -1);
	s.assign(f_num_cells, 1.0);

	p_inv_spacing = 1.0 / (2.2 * p_radius);
	p_num_x = static_cast<int>(std::floor(f_num_x * cell_size * p_inv_spacing)) + 1;
	p_num_y = static_cast<int>(std::floor(f_num_y * cell_size * p_inv_spacing)) + 1;
	int p_num_cells = p_num_x * p_num_y;
	num_cell_particles.assign(p_num_cells, 0);
	first_cell_particle.assign(p_num_cells + 1, 0);

	q_valid = false;
	p_positions.clear();
	p_velocities.clear();
	p_churn.clear();
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

void FLIPWaterSimulation::set_drift_compensation(double p_drift) {
	drift_compensation = std::max(p_drift, 0.0);
}

double FLIPWaterSimulation::get_drift_compensation() const {
	return drift_compensation;
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

void FLIPWaterSimulation::set_particles_per_side(int p_count) {
	particles_per_side = std::clamp(p_count, 1, 4);
	_rebuild_grid();
}

int FLIPWaterSimulation::get_particles_per_side() const {
	return particles_per_side;
}

void FLIPWaterSimulation::set_separation_iterations(int p_iterations) {
	separation_iterations = std::clamp(p_iterations, 0, 16);
}

int FLIPWaterSimulation::get_separation_iterations() const {
	return separation_iterations;
}

void FLIPWaterSimulation::set_drift_pull(bool p_enabled) {
	drift_pull = p_enabled;
}

bool FLIPWaterSimulation::is_drift_pull() const {
	return drift_pull;
}

// Cohesion pulls particles that have drifted apart (between one and two rest spacings) back toward each other,
// a share of the gap per step. Inside a body of water the pulls cancel; at its edge they hold it together, so a
// stream falls as a rope and breaks into blobs instead of fraying into single drops. Viscosity evens each
// particle's velocity toward its neighbours', so water moves as a body. Both 0 leave the sim as it was.
void FLIPWaterSimulation::set_cohesion(double p_amount) {
	cohesion = std::clamp(p_amount, 0.0, 1.0);
}

void FLIPWaterSimulation::set_thin_runoff_ground(bool p_enabled) {
	thin_runoff_ground = p_enabled;
}

bool FLIPWaterSimulation::get_thin_runoff_ground() const {
	return thin_runoff_ground;
}

double FLIPWaterSimulation::get_cohesion() const {
	return cohesion;
}

void FLIPWaterSimulation::set_viscosity(double p_amount) {
	viscosity = std::clamp(p_amount, 0.0, 1.0);
}

double FLIPWaterSimulation::get_viscosity() const {
	return viscosity;
}

// Cells holding less than this fraction of full water's density count as air: their particles fly free under
// gravity (spray, thin sheets, falling streams) instead of being solved as a cell full of water, which made
// sparse water swell into thick blobs. Zero solves every cell that holds a particle.
void FLIPWaterSimulation::set_spray_density(double p_fraction) {
	spray_density = std::clamp(p_fraction, 0.0, 0.9);
}

double FLIPWaterSimulation::get_spray_density() const {
	return spray_density;
}

void FLIPWaterSimulation::set_churn_low(double p_speed) {
	churn_low = std::max(p_speed, 0.0);
}

double FLIPWaterSimulation::get_churn_low() const {
	return churn_low;
}

void FLIPWaterSimulation::set_churn_high(double p_speed) {
	churn_high = std::max(p_speed, 1.0);
}

double FLIPWaterSimulation::get_churn_high() const {
	return churn_high;
}

void FLIPWaterSimulation::set_churn_fade(double p_rate) {
	churn_fade = std::max(p_rate, 0.0);
}

double FLIPWaterSimulation::get_churn_fade() const {
	return churn_fade;
}

void FLIPWaterSimulation::set_wall_drag(double p_rate) {
	wall_drag = std::max(p_rate, 0.0);
}

double FLIPWaterSimulation::get_wall_drag() const {
	return wall_drag;
}

void FLIPWaterSimulation::set_void_fill_cells(int p_cells) {
	void_fill_cells = std::max(p_cells, 0);
}

int FLIPWaterSimulation::get_void_fill_cells() const {
	return void_fill_cells;
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
	double step = cell_size / particles_per_side;
	for (int a = 0; a < particles_per_side; a++) {
		for (int b = 0; b < particles_per_side; b++) {
			_add_particle_internal(Vector2(x + (a + 0.5) * step, y + (b + 0.5) * step), Vector2());
		}
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
	q_valid = false;
	p_positions.clear();
	p_velocities.clear();
	p_churn.clear();
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
	q_valid = false;
	for (size_t i = 0; i < p_positions.size(); i++) {
		p_positions[i] -= offset;
	}
	_remove_particles_outside_domain();
	// The last surface build is kept in local space: move it with the water, or the next build settles toward
	// a field that is off by the shift.
	surf_prev_origin -= offset;
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

void FLIPWaterSimulation::_ensure_query_index() const {
	if (q_valid) {
		return;
	}
	int n = static_cast<int>(p_positions.size());
	auto cell_of = [this](const Vector2 &p) {
		int i = std::clamp(static_cast<int>(std::floor(p.x * f_inv_spacing)), 0, f_num_x - 1);
		int j = std::clamp(static_cast<int>(std::floor(p.y * f_inv_spacing)), 0, f_num_y - 1);
		return i + (f_num_x * j);
	};
	q_first.assign(f_num_cells + 1, 0);
	q_ids.resize(n);
	for (int k = 0; k < n; k++) {
		q_first[cell_of(p_positions[k])]++;
	}
	int sum = 0;
	for (int c = 0; c < f_num_cells; c++) {
		sum += q_first[c];
		q_first[c] = sum;
	}
	q_first[f_num_cells] = sum;
	for (int k = 0; k < n; k++) {
		q_ids[--q_first[cell_of(p_positions[k])]] = k;
	}
	q_valid = true;
}

template <typename F>
void FLIPWaterSimulation::_each_in_rect(const Rect2 &p_rect, F &&p_fn) const {
	if (f_num_cells <= 0) {
		return;
	}
	_ensure_query_index();
	Rect2 r = p_rect.abs();
	int i0 = std::clamp(static_cast<int>(std::floor(r.position.x * f_inv_spacing)), 0, f_num_x - 1);
	int j0 = std::clamp(static_cast<int>(std::floor(r.position.y * f_inv_spacing)), 0, f_num_y - 1);
	int i1 = std::clamp(static_cast<int>(std::floor(r.get_end().x * f_inv_spacing)), 0, f_num_x - 1);
	int j1 = std::clamp(static_cast<int>(std::floor(r.get_end().y * f_inv_spacing)), 0, f_num_y - 1);
	for (int j = j0; j <= j1; j++) {
		for (int i = i0; i <= i1; i++) {
			int c = i + (f_num_x * j);
			for (int k = q_first[c]; k < q_first[c + 1]; k++) {
				p_fn(q_ids[k]);
			}
		}
	}
}

template <typename F>
void FLIPWaterSimulation::_each_in_circle(const Vector2 &p_center, double p_radius, F &&p_fn) const {
	double r2 = p_radius * p_radius;
	_each_in_rect(Rect2(p_center - Vector2(p_radius, p_radius), Vector2(p_radius, p_radius) * 2.0), [&](int id) {
		if (p_positions[id].distance_squared_to(p_center) <= r2) {
			p_fn(id);
		}
	});
}

bool FLIPWaterSimulation::_add_particle_internal(const Vector2 &pos, const Vector2 &vel) {
	if (static_cast<int>(p_positions.size()) >= max_particles) {
		return false;
	}
	q_valid = false;
	p_positions.push_back(pos);
	p_velocities.push_back(vel);
	p_churn.push_back(0.0f);
	return true;
}

void FLIPWaterSimulation::_remove_particles_outside_domain() {
	q_valid = false;
	Rect2 domain = get_domain_rect();
	size_t kept = 0;
	for (size_t i = 0; i < p_positions.size(); i++) {
		if (domain.has_point(p_positions[i])) {
			p_positions[kept] = p_positions[i];
			p_velocities[kept] = p_velocities[i];
			p_churn[kept] = p_churn[i];
			kept++;
		}
	}
	p_positions.resize(kept);
	p_velocities.resize(kept);
	p_churn.resize(kept);
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
	_each_in_circle(p_center, p_radius, [&](int i) { p_velocities[i] += p_velocity; });
}

int FLIPWaterSimulation::remove_particles_in_circle(const Vector2 &p_center, double p_radius) {
	// Found through the index; only an actual removal pays for compacting the arrays.
	int found = 0;
	_each_in_circle(p_center, p_radius, [&](int) { found++; });
	if (found == 0) {
		return 0;
	}
	q_valid = false;
	double r2 = p_radius * p_radius;
	size_t kept = 0;
	for (size_t i = 0; i < p_positions.size(); i++) {
		if (p_positions[i].distance_squared_to(p_center) > r2) {
			p_positions[kept] = p_positions[i];
			p_velocities[kept] = p_velocities[i];
			p_churn[kept] = p_churn[i];
			kept++;
		}
	}
	int removed = static_cast<int>(p_positions.size() - kept);
	p_positions.resize(kept);
	p_velocities.resize(kept);
	p_churn.resize(kept);
	return removed;
}

int FLIPWaterSimulation::count_particles_in_circle(const Vector2 &p_center, double p_radius) const {
	int count = 0;
	_each_in_circle(p_center, p_radius, [&](int) { count++; });
	return count;
}

int FLIPWaterSimulation::count_particles_in_rect(const Rect2 &p_rect) const {
	Rect2 rect = p_rect.abs();
	int count = 0;
	_each_in_rect(rect, [&](int i) {
		if (rect.has_point(p_positions[i])) {
			count++;
		}
	});
	return count;
}

Vector2 FLIPWaterSimulation::get_average_velocity_in_circle(const Vector2 &p_center, double p_radius) const {
	Vector2 sum;
	int count = 0;
	_each_in_circle(p_center, p_radius, [&](int i) {
		sum += p_velocities[i];
		count++;
	});
	return count > 0 ? sum / count : Vector2();
}

void FLIPWaterSimulation::apply_radial_impulse(const Vector2 &p_center, double p_radius, double p_strength) {
	_each_in_circle(p_center, p_radius, [&](int i) {
		Vector2 offset = p_positions[i] - p_center;
		double d2 = offset.length_squared();
		if (d2 > 0.0) {
			double d = std::sqrt(d2);
			p_velocities[i] += offset / d * (p_strength * (1.0 - d / p_radius));
		}
	});
}

// Water that enters or leaves the grid inside a circle this frame, `rate` px^2 of it per second (negative drains):
// the pressure solve is told the fluid cells there gain or lose that much volume, so the water round them flows in
// to fill what is drained, or out to make room for what is added, in the same step. Taking particles out (or adding
// them) without it left a hole the solver did not know to fill, or a crowd it threw apart. The caller still removes
// or adds the particles themselves, at the same rate. Lasts until the end of the next step().
void FLIPWaterSimulation::add_volume_source(const Vector2 &p_center, double p_radius, double p_rate) {
	if (p_rate == 0.0 || p_radius <= 0.0) {
		return;
	}
	volume_sources.push_back({ p_center, p_radius, p_rate });
}

void FLIPWaterSimulation::_apply_volume_sources() {
	for (const VolumeSource &src : volume_sources) {
		int i0 = std::clamp(static_cast<int>(std::floor((src.center.x - src.radius) * f_inv_spacing)), 1, f_num_x - 2);
		int i1 = std::clamp(static_cast<int>(std::floor((src.center.x + src.radius) * f_inv_spacing)), 1, f_num_x - 2);
		int j0 = std::clamp(static_cast<int>(std::floor((src.center.y - src.radius) * f_inv_spacing)), 1, f_num_y - 2);
		int j1 = std::clamp(static_cast<int>(std::floor((src.center.y + src.radius) * f_inv_spacing)), 1, f_num_y - 2);
		double r2 = src.radius * src.radius;
		source_cells.clear();
		for (int j = j0; j <= j1; j++) {
			for (int i = i0; i <= i1; i++) {
				int idx = i + (j * f_num_x);
				if (cell_types[idx] != FLUID_CELL) {
					continue;
				}
				Vector2 c((i + 0.5) * cell_size, (j + 0.5) * cell_size);
				if (c.distance_squared_to(src.center) <= r2) {
					source_cells.push_back(idx);
				}
			}
		}
		if (source_cells.empty()) {
			// A circle smaller than a cell: the fluid cell it sits in, if any.
			int i = std::clamp(static_cast<int>(std::floor(src.center.x * f_inv_spacing)), 1, f_num_x - 2);
			int j = std::clamp(static_cast<int>(std::floor(src.center.y * f_inv_spacing)), 1, f_num_y - 2);
			int idx = i + (j * f_num_x);
			if (cell_types[idx] == FLUID_CELL) {
				source_cells.push_back(idx);
			}
		}
		if (source_cells.empty()) {
			continue;
		}
		// The solve's divergence is a cell's net outflow in px/s across its faces, so a volume rate per cell divides
		// by the cell's size once.
		double per_cell = src.rate / (static_cast<double>(source_cells.size()) * cell_size);
		for (int idx : source_cells) {
			drift_source[idx] += per_cell;
		}
	}
}

// Writes one 2D instance transform per particle (optionally only those moving at least `min_speed`),
// so custom particle rendering takes a single draw call and no per-particle script work.
// `every` > 1 keeps only every n-th qualifying particle, for sparse markers such as bubbles. `min_churn`
// keeps only particles at least that churned. A multimesh that uses colours gets each particle's churn in
// its colour (r, g and b), for drawing white water.
// Surface reconstruction (Zhu & Bridson 2005): each sample's distance to the kernel-weighted average of the particles
// near it, less a radius. Flat-packed particles give a flat surface, not a row of blobs. Particles touching a solid face
// are mirrored across it, so the field meets walls square (no meniscus) and fills the gap particles keep from walls.
// That distance is only right near the surface (it bottoms out at -radius), so a distance transform then gives depth
// below it, through solids as if they were water, so walls don't count as surface.
// Output is RGBA8, one texel per sample: red the signed distance over +-band, green the weighted churn, blue the depth
// below the open surface over 0..depth px (in air: how closed the sample's pocket is), alpha the weighted particle speed over 0..max_speed.
PackedByteArray FLIPWaterSimulation::build_surface(const Vector2 &p_origin, const Vector2i &p_samples, double p_spacing, double p_kernel, double p_radius, double p_band, double p_ghost_reach, int p_smooth, double p_level, double p_level_tol, double p_settle, double p_depth, double p_pocket) {
	int w = std::max(p_samples.x, 1);
	int ht = std::max(p_samples.y, 1);
	size_t n = static_cast<size_t>(w) * ht;
	surf_w.assign(n, 0.0f);
	surf_x.assign(n, 0.0f);
	surf_y.assign(n, 0.0f);
	surf_c.assign(n, 0.0f);
	surf_s.assign(n, 0.0f);
	surf_phi.resize(n);
	surf_tmp.resize(n);
	surf_solid.resize(n);
	surf_nx = w;
	surf_ny = ht;
	surf_sp = static_cast<float>(std::max(p_spacing, 0.5));
	surf_kernel = static_cast<float>(std::max(p_kernel, p_spacing));
	surf_reach = static_cast<float>(p_ghost_reach);
	surf_origin = p_origin;
	float sp = surf_sp;


	// The particles in reach, with the solid faces they touch.
	double h = cell_size;
	auto solid = [&](int ci, int cj) {
		return !_is_cell_in_bounds(ci, cj) || s[ci + f_num_x * cj] == 0.0;
	};
	float margin = surf_kernel + surf_reach + 2.0f * sp;
	float lo_x = p_origin.x - margin, hi_x = p_origin.x + w * sp + margin;
	float lo_y = p_origin.y - margin, hi_y = p_origin.y + ht * sp + margin;
	const float none = std::numeric_limits<float>::quiet_NaN();
	const float ceiling_reach = static_cast<float>(1.25 * this->p_radius);
	const double loaded_pressure = std::max(DENSITY * static_cast<double>(gravity.length()) * cell_size, 1e-6);
	surf_parts.clear();
	for (size_t p = 0; p < p_positions.size(); p++) {
		float px = p_positions[p].x;
		float py = p_positions[p].y;
		if (px < lo_x || px > hi_x || py < lo_y || py > hi_y) {
			continue;
		}
		SurfParticle sp_ = { px, py, p_churn[p], static_cast<float>(p_velocities[p].length()), none, none };
		if (surf_reach > 0.0f) {
			int ci = static_cast<int>(std::floor(px / h));
			int cj = static_cast<int>(std::floor(py / h));
			if (px - ci * h < surf_reach && solid(ci - 1, cj)) {
				sp_.face_x = static_cast<float>(ci * h);
			} else if ((ci + 1) * h - px < surf_reach && solid(ci + 1, cj)) {
				sp_.face_x = static_cast<float>((ci + 1) * h);
			}
			// A face overhead mirrors water according to the pressure under it. Water under load (a cell or more of
			// water standing above it somewhere) is pushed against the face, and the gap the particles leave there
			// is not air, so it is mirrored across a whole cell. Water under no load hangs below the block with a
			// surface of its own: mirroring it would draw it reaching up to hug the block, so there only particles
			// pressed against the face are mirrored.
			float overhead = surf_reach;
			if (solid(ci, cj - 1) && _is_cell_in_bounds(ci, cj)) {
				float loaded = static_cast<float>(std::clamp(pressure[ci + f_num_x * cj] / loaded_pressure - 0.5, 0.0, 1.0));
				overhead = ceiling_reach + (std::max(surf_reach, static_cast<float>(h)) - ceiling_reach) * loaded;
			}
			if (py - cj * h < overhead && solid(ci, cj - 1)) {
				sp_.face_y = static_cast<float>(cj * h);
			} else if ((cj + 1) * h - py < surf_reach && solid(ci, cj + 1)) {
				sp_.face_y = static_cast<float>((cj + 1) * h);
			}
		}
		surf_parts.push_back(sp_);
	}

	// Rows are independent: each chunk splats every particle near its rows into those rows only.
	_run_parallel(ht, 8, &FLIPWaterSimulation::_surface_rows);

	float *W = surf_w.data(), *C = surf_c.data();
	surf_band = static_cast<float>(std::max(p_band, 0.5));
	surf_radius = static_cast<float>(p_radius);
	float band = surf_band;
	_run_parallel(ht, 8, &FLIPWaterSimulation::_surface_phi_rows);

	// [1 2 1] passes, separable: evens out the particles' small jitter along the surface.
	for (int pass = 0; pass < p_smooth; pass++) {
		_run_parallel(ht, 8, &FLIPWaterSimulation::_surface_smooth_x_rows);
		_run_parallel(ht, 8, &FLIPWaterSimulation::_surface_smooth_y_rows);
	}

	// Pockets: air samples (outside water and solids) are grouped into connected patches. A patch that reaches the
	// grid's edge is open air. An enclosed one is a pocket, and how open it counts grows with its area, from a
	// quarter of `pocket` px^2 to all of it. Small pockets (a gap in a corner, a vortex core just forming) are drawn
	// filled in, and they open up as they grow rather than popping in with a full rim and glow.
	std::vector<float> &d = surf_tmp;
	const float depth_range = std::max(static_cast<float>(p_depth), band);
	const float far = depth_range + 2.0f * sp;
	float pocket_full = std::max(static_cast<float>(p_pocket), 1.0f) / (sp * sp);
	float pocket_none = 0.25f * pocket_full;
	std::vector<int32_t> &stack = surf_stack;
	std::vector<int32_t> &patch = surf_patch;
	std::vector<uint8_t> &seen = surf_seen;
	seen.assign(n, 0);
	for (size_t idx = 0; idx < n; idx++) {
		bool air = surf_phi[idx] >= 0.0f && !surf_solid[idx];
		d[idx] = far;
		if (!air) {
			seen[idx] = 1;
		}
	}
	for (size_t start = 0; start < n; start++) {
		if (seen[start]) {
			continue;
		}
		bool open = false;
		patch.clear();
		stack.clear();
		stack.push_back(static_cast<int32_t>(start));
		seen[start] = 1;
		while (!stack.empty()) {
			int32_t idx = stack.back();
			stack.pop_back();
			patch.push_back(idx);
			int i = idx % w;
			int j = idx / w;
			if (i == 0 || j == 0 || i == w - 1 || j == ht - 1) {
				open = true;
			}
			if (i > 0 && !seen[idx - 1]) {
				seen[idx - 1] = 1;
				stack.push_back(idx - 1);
			}
			if (i < w - 1 && !seen[idx + 1]) {
				seen[idx + 1] = 1;
				stack.push_back(idx + 1);
			}
			if (j > 0 && !seen[idx - w]) {
				seen[idx - w] = 1;
				stack.push_back(idx - w);
			}
			if (j < ht - 1 && !seen[idx + w]) {
				seen[idx + w] = 1;
				stack.push_back(idx + w);
			}
		}
		float openness = 1.0f;
		if (!open) {
			float t = std::clamp((patch.size() - pocket_none) / std::max(pocket_full - pocket_none, 1.0f), 0.0f, 1.0f);
			openness = t * t * (3.0f - 2.0f * t);
		}
		float cost = (1.0f - openness) * depth_range;
		for (int32_t idx : patch) {
			d[idx] = cost;
		}
	}

	// Depth: chamfer distance from air into water and solids, on top of what the air's patch costs. Sideways steps
	// cost less than vertical ones, so depth mostly means "how far below the surface": the shading runs in level
	// bands that reach in under overhangs, instead of fanning out from wherever the surface meets a block.
	const float side = 0.45f;
	float dx1 = sp * side;
	float dy1 = sp;
	float dd = sp * std::sqrt(1.0f + side * side);
	for (int j = 1; j < ht; j++) {
		float *row = d.data() + static_cast<size_t>(j) * w;
		const float *up = row - w;
		for (int i = 1; i < w - 1; i++) {
			float v = std::min(std::min(row[i - 1] + dx1, up[i] + dy1), std::min(up[i - 1] + dd, up[i + 1] + dd));
			row[i] = std::min(row[i], v);
		}
	}
	for (int j = ht - 2; j >= 0; j--) {
		float *row = d.data() + static_cast<size_t>(j) * w;
		const float *down = row + w;
		for (int i = w - 2; i >= 1; i--) {
			float v = std::min(std::min(row[i + 1] + dx1, down[i] + dy1), std::min(down[i + 1] + dd, down[i - 1] + dd));
			row[i] = std::min(row[i], v);
		}
	}
	for (size_t idx = 0; idx < n; idx++) {
		if (surf_phi[idx] < 0.0f) {
			// Near the surface the reconstruction is exact; deeper, the transform's depth (half a sample short, so
			// the two meet without a step). The signed distance only ever needs true distance, so pockets and the
			// cheaper sideways steps are left out of it by never letting the transform shorten it.
			float depth = std::max(-surf_phi[idx], d[idx] - 0.5f * sp);
			surf_phi[idx] = -std::min(std::max(-surf_phi[idx], std::min(depth, band)), band);
			d[idx] = std::min(depth, depth_range);
		} else {
			// In air this holds how closed the sample's pocket is, in the same px: zero in open air.
			d[idx] = std::min(d[idx], depth_range);
		}
	}

	// The depth goes out in blue, blurred a little so it shades smoothly, and eased toward last frame's so shading
	// and pockets fade in and out over a few frames. Air that was water a frame ago takes its value at once:
	// easing there would leave a ghost of the water behind a surface that is falling away.
	surf_depth.resize(n);
	bool have_prev = surf_prev_sp == sp && !surf_prev.empty() && surf_depth_prev.size() == surf_prev.size();
	int shift_i = static_cast<int>(std::lround((p_origin.x - surf_prev_origin.x) / sp));
	int shift_j = static_cast<int>(std::lround((p_origin.y - surf_prev_origin.y) / sp));
	const float ease = 0.3f;
	for (int j = 0; j < ht; j++) {
		for (int i = 0; i < w; i++) {
			size_t idx = static_cast<size_t>(j) * w + i;
			float sum = 0.0f;
			for (int dj = -1; dj <= 1; dj++) {
				int jj = std::clamp(j + dj, 0, ht - 1);
				for (int di = -1; di <= 1; di++) {
					sum += d[static_cast<size_t>(jj) * w + std::clamp(i + di, 0, w - 1)];
				}
			}
			float now = std::clamp(sum / (9.0f * depth_range), 0.0f, 1.0f);
			int pi = i + shift_i;
			int pj = j + shift_j;
			if (have_prev && pi >= 0 && pi < surf_prev_nx && pj >= 0 && pj < surf_prev_ny) {
				size_t pidx = static_cast<size_t>(pj) * surf_prev_nx + pi;
				bool air = surf_phi[idx] >= 0.0f && !surf_solid[idx];
				if (!(air && surf_prev[pidx] < 0.0f)) {
					float before = surf_depth_prev[pidx];
					now = before + (now - before) * ease;
				}
			}
			surf_depth[idx] = now;
		}
	}
	surf_depth_prev = surf_depth;

	// Levelling: each sample averages its row out to `level` px either side, stopping at solids, over neighbours whose
	// distance is within `level_tol` of what the row's local slope predicts. Along a surface, neighbours differ from
	// that by its dips, bumps and particle-row steps, so those even out into a smooth line while broad waves stay.
	// Across a drop's edge they differ by more and are left out, so splashes keep their shape.
	// Both passes are for calm water only: an upward-facing surface with water at least ~16 px deep under it. Streams,
	// drops and sheets are left exactly as the particles put them; levelling them across frames whittles them into
	// shards.
	{
		surf_flat.assign(n, 0.0f);
		int below = std::max(1, static_cast<int>(std::lround(16.0f / sp)));
		// The slope is read over ±3 samples, so a particle-row step (which tilts the surface right at it) still counts as flat.
		int s3 = 3;
		for (int j = s3; j < ht - s3; j++) {
			for (int i = s3; i < w - s3; i++) {
				size_t idx = static_cast<size_t>(j) * w + i;
				if (surf_solid[idx] || j + below >= ht) {
					continue;
				}
				float gx = (surf_phi[idx + s3] - surf_phi[idx - s3]) / (2.0f * s3);
				float gy = (surf_phi[idx + static_cast<size_t>(s3) * w] - surf_phi[idx - static_cast<size_t>(s3) * w]) / (2.0f * s3);
				float g = std::sqrt(gx * gx + gy * gy);
				if (g < 0.2f * sp) {
					continue;
				}
				float up = std::clamp((-gy / g - 0.5f) / 0.3f, 0.0f, 1.0f);
				float deep = std::clamp((-surf_phi[idx + static_cast<size_t>(below) * w] - 4.0f) / 8.0f, 0.0f, 1.0f);
				surf_flat[idx] = up * up * (3.0f - 2.0f * up) * deep * deep * (3.0f - 2.0f * deep);
			}
		}
	}
	surf_level = static_cast<int>(std::lround(std::max(p_level, 0.0) / sp));
	surf_level_tol = static_cast<float>(std::max(p_level_tol, 0.01));
	if (surf_level > 0) {
		_run_parallel(ht, 8, &FLIPWaterSimulation::_surface_level_rows);
		surf_phi.swap(surf_tmp);
	}

	// Settling: blends toward last frame's field where the two are close, which calms the particles' frame-to-frame
	// jitter. Where they differ by more than the tolerance (water that moved), this frame wins outright, so nothing trails.
	float settle = static_cast<float>(std::clamp(p_settle, 0.0, 0.95));
	if (settle > 0.0f && surf_prev_sp == sp && !surf_prev.empty()) {
		int di = static_cast<int>(std::lround((p_origin.x - surf_prev_origin.x) / sp));
		int dj = static_cast<int>(std::lround((p_origin.y - surf_prev_origin.y) / sp));
		float tol = surf_level_tol;
		for (int j = 0; j < ht; j++) {
			int pj = j + dj;
			if (pj < 0 || pj >= surf_prev_ny) {
				continue;
			}
			for (int i = 0; i < w; i++) {
				int pi = i + di;
				if (pi < 0 || pi >= surf_prev_nx) {
					continue;
				}
				size_t idx = static_cast<size_t>(j) * w + i;
				float prev = surf_prev[static_cast<size_t>(pj) * surf_prev_nx + pi];
				float diff = (prev - surf_phi[idx]) / tol;
				float keep = settle * surf_flat[idx] * std::max(0.0f, 1.0f - diff * diff);
				surf_phi[idx] += (prev - surf_phi[idx]) * keep;
			}
		}
	}
	surf_prev = surf_phi;
	surf_prev_origin = p_origin;
	surf_prev_nx = w;
	surf_prev_ny = ht;
	surf_prev_sp = sp;

	float speed_range = max_speed > 0.0 ? static_cast<float>(max_speed) : 1000.0f;
	PackedByteArray out;
	out.resize(static_cast<int64_t>(n) * 4);
	uint8_t *o = out.ptrw();
	for (size_t idx = 0; idx < n; idx++) {
		float red = 0.5f + surf_phi[idx] / (2.0f * band);
		float churn = W[idx] > 1e-6f ? C[idx] : 0.0f;
		o[idx * 4 + 0] = static_cast<uint8_t>(std::lround(std::clamp(red, 0.0f, 1.0f) * 255.0f));
		o[idx * 4 + 1] = static_cast<uint8_t>(std::lround(std::clamp(churn, 0.0f, 1.0f) * 255.0f));
		float speed = W[idx] > 1e-6f ? surf_s[idx] / speed_range : 0.0f;
		o[idx * 4 + 2] = static_cast<uint8_t>(std::lround(surf_depth[idx] * 255.0f));
		o[idx * 4 + 3] = static_cast<uint8_t>(std::lround(std::clamp(speed, 0.0f, 1.0f) * 255.0f));
	}
	return out;
}

void FLIPWaterSimulation::_surface_phi_rows(int begin, int end) {
	int w = surf_nx;
	float sp = surf_sp;
	float ox = surf_origin.x;
	float oy = surf_origin.y;
	float band = surf_band;
	float r = surf_radius;
	double inv_h = 1.0 / cell_size;
	float *W = surf_w.data(), *X = surf_x.data(), *Y = surf_y.data(), *C = surf_c.data();
	for (int j = begin; j < end; j++) {
		float sy = oy + (j + 0.5f) * sp;
		int cj = static_cast<int>(std::floor(sy * inv_h));
		for (int i = 0; i < w; i++) {
			size_t idx = static_cast<size_t>(j) * w + i;
			float sx = ox + (i + 0.5f) * sp;
			float phi = band;
			if (W[idx] > 1e-6f) {
				float cx = X[idx] / W[idx] - sx;
				float cy = Y[idx] / W[idx] - sy;
				phi = std::sqrt(cx * cx + cy * cy) - r;
				C[idx] /= W[idx];
				surf_s[idx] /= W[idx];
			}
			surf_phi[idx] = std::clamp(phi, -band, band);
			int ci = static_cast<int>(std::floor(sx * inv_h));
			surf_solid[idx] = (!_is_cell_in_bounds(ci, cj) || s[ci + f_num_x * cj] == 0.0) ? 1 : 0;
		}
	}
}

void FLIPWaterSimulation::_surface_smooth_x_rows(int begin, int end) {
	int w = surf_nx;
	for (int j = begin; j < end; j++) {
		const float *src = surf_phi.data() + static_cast<size_t>(j) * w;
		float *dst = surf_tmp.data() + static_cast<size_t>(j) * w;
		for (int i = 0; i < w; i++) {
			dst[i] = 0.25f * (src[i > 0 ? i - 1 : i] + 2.0f * src[i] + src[i < w - 1 ? i + 1 : i]);
		}
	}
}

void FLIPWaterSimulation::_surface_smooth_y_rows(int begin, int end) {
	int w = surf_nx;
	int ht = surf_ny;
	for (int j = begin; j < end; j++) {
		const float *mid = surf_tmp.data() + static_cast<size_t>(j) * w;
		const float *up = surf_tmp.data() + static_cast<size_t>(j > 0 ? j - 1 : j) * w;
		const float *down = surf_tmp.data() + static_cast<size_t>(j < ht - 1 ? j + 1 : j) * w;
		float *dst = surf_phi.data() + static_cast<size_t>(j) * w;
		for (int i = 0; i < w; i++) {
			dst[i] = 0.25f * (up[i] + 2.0f * mid[i] + down[i]);
		}
	}
}

void FLIPWaterSimulation::_surface_level_rows(int begin, int end) {
	int w = surf_nx;
	int k = surf_level;
	float inv_tol = 1.0f / surf_level_tol;
	float inv_k = 1.0f / (k + 1);
	for (int j = begin; j < end; j++) {
		const float *src = surf_phi.data() + static_cast<size_t>(j) * w;
		const uint8_t *sol = surf_solid.data() + static_cast<size_t>(j) * w;
		float *dst = surf_tmp.data() + static_cast<size_t>(j) * w;
		const float *flat = surf_flat.data() + static_cast<size_t>(j) * w;
		for (int i = 0; i < w; i++) {
			float self = src[i];
			if (sol[i]) {
				dst[i] = self;
				continue;
			}
			// The row's slope across the whole reach (wide, so a step in the middle barely tilts it), so a sloping surface
			// compares neighbours with where the slope puts them.
			int l = i, r = i;
			for (int o = 1; o <= k && l > 0 && !sol[l - 1]; o++) {
				l--;
			}
			for (int o = 1; o <= k && r < w - 1 && !sol[r + 1]; o++) {
				r++;
			}
			float slope = r > l ? (src[r] - src[l]) / (r - l) : 0.0f;
			// Averages how far neighbours sit off that slope, so a one-sided average (by a wall) doesn't drag it.
			float sum = 0.0f;
			float wsum = 1.0f;
			for (int dir = -1; dir <= 1; dir += 2) {
				for (int o = 1; o <= k; o++) {
					int ni = i + dir * o;
					if (ni < 0 || ni >= w || sol[ni]) {
						break;
					}
					float s_ = 1.0f - (o * inv_k) * (o * inv_k);
					float off = src[ni] - self - slope * (dir * o);
					float r_ = off * inv_tol;
					float r2 = 1.0f - r_ * r_;
					if (r2 <= 0.0f) {
						continue;
					}
					float wt = s_ * s_ * r2 * r2;
					sum += off * wt;
					wsum += wt;
				}
			}
			dst[i] = self + flat[i] * sum / wsum;
		}
	}
}

void FLIPWaterSimulation::_surface_rows(int begin, int end) {
	int w = surf_nx;
	float sp = surf_sp;
	float inv = 1.0f / sp;
	float R = surf_kernel;
	float R2 = R * R;
	float invR2 = 1.0f / R2;
	float ox = surf_origin.x;
	float oy = surf_origin.y;
	float lo_x = ox - R, hi_x = ox + w * sp + R;
	float lo_y = oy + begin * sp - R, hi_y = oy + end * sp + R;
	float *W = surf_w.data(), *X = surf_x.data(), *Y = surf_y.data(), *C = surf_c.data();

	float *S = surf_s.data();
	double h = cell_size;
	auto solid = [&](int ci, int cj) {
		return !_is_cell_in_bounds(ci, cj) || s[ci + f_num_x * cj] == 0.0;
	};
	// A ghost (a particle mirrored across a solid face) stands for water pressed against that face, so it only counts
	// where the face is: at samples whose row (a side face) or column (a floor or ceiling) runs into the block it was
	// mirrored across. Anywhere else it would pull the surface toward a block the water isn't touching: up the sides
	// of a pillar the water passes under, or over the lip of a ledge. `wall_x`/`wall_y` are the cell column/row of
	// that block (-1 none).
	auto splat = [&](float px, float py, float churn, float speed, int wall_x = -1, int wall_y = -1) {
		if (px < lo_x || px > hi_x || py < lo_y || py > hi_y) {
			return;
		}
		int i0 = std::max(0, static_cast<int>(std::ceil((px - R - ox) * inv - 0.5f)));
		int i1 = std::min(w - 1, static_cast<int>(std::floor((px + R - ox) * inv - 0.5f)));
		int j0 = std::max(begin, static_cast<int>(std::ceil((py - R - oy) * inv - 0.5f)));
		int j1 = std::min(end - 1, static_cast<int>(std::floor((py + R - oy) * inv - 0.5f)));
		for (int j = j0; j <= j1; j++) {
			float dy = oy + (j + 0.5f) * sp - py;
			float dy2 = dy * dy;
			size_t row = static_cast<size_t>(j) * w;
			int cell_j = static_cast<int>(std::floor((oy + (j + 0.5f) * sp) / h));
			if (wall_x >= 0 && !solid(wall_x, cell_j)) {
				continue;
			}
			for (int i = i0; i <= i1; i++) {
				float dx = ox + (i + 0.5f) * sp - px;
				float d2 = dx * dx + dy2;
				if (d2 >= R2) {
					continue;
				}
				if (wall_y >= 0 && !solid(static_cast<int>(std::floor((ox + (i + 0.5f) * sp) / h)), wall_y)) {
					continue;
				}
				float q = 1.0f - d2 * invR2;
				float k = q * q * q;
				size_t idx = row + i;
				W[idx] += k;
				X[idx] += k * px;
				Y[idx] += k * py;
				C[idx] += k * churn;
				S[idx] += k * speed;
			}
		}
	};
	for (const SurfParticle &q : surf_parts) {
		splat(q.x, q.y, q.churn, q.speed);
		bool hx = !std::isnan(q.face_x);
		bool hy = !std::isnan(q.face_y);
		float gx = hx ? 2.0f * q.face_x - q.x : q.x;
		float gy = hy ? 2.0f * q.face_y - q.y : q.y;
		// The block each face belongs to: on the ghost's side of it.
		int wx = hx ? static_cast<int>(std::floor(q.face_x / h)) - (q.x > q.face_x ? 1 : 0) : -1;
		int wy = hy ? static_cast<int>(std::floor(q.face_y / h)) - (q.y > q.face_y ? 1 : 0) : -1;
		if (hx) {
			splat(gx, q.y, q.churn, q.speed, wx, -1);
		}
		if (hy) {
			splat(q.x, gy, q.churn, q.speed, -1, wy);
		}
		if (hx && hy) {
			splat(gx, gy, q.churn, q.speed, wx, wy);
		}
	}
}

// Up to `probes` random particles in the rect, each kept with probability equal to its churn (if at least min_churn).
PackedVector2Array FLIPWaterSimulation::sample_churned(const Rect2 &p_rect, int p_probes, double p_min_churn) {
	PackedVector2Array out;
	size_t count = p_positions.size();
	if (count == 0) {
		return out;
	}
	std::uniform_real_distribution<float> unit(0.0f, 1.0f);
	for (int k = 0; k < p_probes; k++) {
		size_t i = static_cast<size_t>(rng()) % count;
		float churn = p_churn[i];
		if (churn >= p_min_churn && churn > 0.0f && unit(rng) < churn && p_rect.has_point(p_positions[i])) {
			out.push_back(p_positions[i]);
		}
	}
	return out;
}

// The grid's velocity at a point, interpolated as particles read it.
Vector2 FLIPWaterSimulation::sample_velocity(const Vector2 &p_pos) const {
	if (f_num_x < 2 || f_num_y < 2) {
		return Vector2();
	}
	double h = cell_size;
	double h1 = f_inv_spacing;
	double h2 = 0.5 * h;
	Vector2 result;
	for (int component = 0; component < 2; component++) {
		double dx = component == 0 ? 0.0 : h2;
		double dy = component == 0 ? h2 : 0.0;
		const std::vector<double> &f = component == 0 ? u : v;
		double x = std::clamp(static_cast<double>(p_pos.x), h, f_num_x * h);
		double y = std::clamp(static_cast<double>(p_pos.y), h, f_num_y * h);
		int x0 = std::min(static_cast<int>(std::floor((x - dx) * h1)), f_num_x - 1);
		double tx = ((x - dx) - x0 * h) * h1;
		int x1 = std::min(x0 + 1, f_num_x - 1);
		int y0 = std::min(static_cast<int>(std::floor((y - dy) * h1)), f_num_y - 1);
		double ty = ((y - dy) - y0 * h) * h1;
		int y1 = std::min(y0 + 1, f_num_y - 1);
		double sx = 1.0 - tx;
		double sy = 1.0 - ty;
		double wts[4] = { sx * sy, tx * sy, tx * ty, sx * ty };
		int nr[4] = { x0 + (f_num_x * y0), x1 + (f_num_x * y0), x1 + (f_num_x * y1), x0 + (f_num_x * y1) };
		int offset = component == 1 ? f_num_x : 1;
		double total = 0.0;
		double sum = 0.0;
		for (int k = 0; k < 4; k++) {
			if (nr[k] - offset >= 0 && cell_types[nr[k]] != SOLID_CELL && cell_types[nr[k] - offset] != SOLID_CELL) {
				total += wts[k];
				sum += wts[k] * f[nr[k]];
			}
		}
		double value = total > 0.0 ? sum / total : 0.0;
		if (component == 0) {
			result.x = value;
		} else {
			result.y = value;
		}
		if (total > 0.0) {
			// As the particles that read it are given (see _grid_to_particles_range): still water reads zero.
			result += Vector2(component == 0 ? gravity.x : 0.0, component == 1 ? gravity.y : 0.0) * step_dt;
		}
	}
	return result;
}

// `mirror` > 0 also writes a copy of each particle within that many px of a solid face, mirrored across it (and one
// mirrored across both faces in a corner), as build_surface does with its ghosts: water drawn by summing splats then
// meets a wall square instead of stopping short of it where particles keep their gap and have neighbours on one side
// only. A ceiling mirrors only water pressed up against it (under load, or filling the space under it), so water
// hanging below a block keeps a surface of its own instead of reaching up to hug it.
//
// `ground` > 0, on a mesh that uses colours, tells the mesh's shader how each particle lies, so that it can draw
// loose water by its shape and not as a ball. Green is 1 for a particle within `ground` px of a floor under it, and
// blue is then that distance over `ground`; its quad is left square to the world, unturned and unstretched. Green is
// 0 for any other particle, and blue is then its speed over 600 px/s. Green is a half for the copies and the
// fillers, which are there for the summed density only. Red stays the churn.
void FLIPWaterSimulation::update_multimesh(const Ref<MultiMesh> &p_multimesh, double p_scale, double p_min_speed, double p_stretch, double p_max_stretch, int p_every, double p_min_churn, double p_mirror, double p_ground) const {
	ERR_FAIL_COND(p_multimesh.is_null());
	ERR_FAIL_COND_MSG(p_multimesh->get_transform_format() != MultiMesh::TRANSFORM_2D, "MultiMesh must use TRANSFORM_2D.");
	ERR_FAIL_COND_MSG(p_multimesh->is_using_custom_data(), "MultiMesh must not use custom data.");
	bool colors = p_multimesh->is_using_colors();
	int stride = colors ? 12 : 8;
	float min_churn = static_cast<float>(p_min_churn);
	double min2 = p_min_speed * p_min_speed;
	size_t every = static_cast<size_t>(std::max(p_every, 1));

	struct Instance {
		Vector2 pos, vel;
		float churn;
		float lies = 0.5f; // the colour's green: 1 on a floor, 0 in the air, a half for a copy
		float floor = 0.0f; // px down to the floor, for one that lies on it
	};
	thread_local std::vector<Instance> instances;
	instances.clear();
	double h = cell_size;
	float reach = static_cast<float>(p_mirror);
	auto solid = [&](int ci, int cj) {
		return !_is_cell_in_bounds(ci, cj) || s[ci + f_num_x * cj] == 0.0;
	};
	const float ceiling_reach = static_cast<float>(1.25 * p_radius);
	const double loaded_pressure = std::max(DENSITY * static_cast<double>(gravity.length()) * cell_size, 1e-6);
	for (size_t i = 0; i < p_positions.size(); i += every) {
		if (p_velocities[i].length_squared() < min2 || p_churn[i] < min_churn) {
			continue;
		}
		Vector2 p = p_positions[i];
		Vector2 v = p_velocities[i];
		int ci = static_cast<int>(std::floor(p.x / h));
		int cj = static_cast<int>(std::floor(p.y / h));
		float to_floor = static_cast<float>((cj + 1) * h - p.y);
		bool on_floor = p_ground > 0.0 && to_floor < p_ground && solid(ci, cj + 1) && !solid(ci, cj);
		instances.push_back({ p, v, p_churn[i], on_floor ? 1.0f : 0.0f, on_floor ? to_floor : 0.0f });
		if (reach <= 0.0f) {
			continue;
		}
		double fx = std::numeric_limits<double>::quiet_NaN();
		double fy = fx;
		if (p.x - ci * h < reach && solid(ci - 1, cj)) {
			fx = ci * h;
		} else if ((ci + 1) * h - p.x < reach && solid(ci + 1, cj)) {
			fx = (ci + 1) * h;
		}
		float overhead = reach;
		if (solid(ci, cj - 1) && _is_cell_in_bounds(ci, cj)) {
			// Water is pressed up against a ceiling when it is under load, or when there is water under it: water held
			// up in a sealed cave has no load. A drip hanging under a block has air below it, and keeps its own surface.
			float loaded = static_cast<float>(std::clamp(pressure[ci + f_num_x * cj] / loaded_pressure - 0.5, 0.0, 1.0));
			if (_is_cell_in_bounds(ci, cj + 1) && cell_types[ci + f_num_x * (cj + 1)] == FLUID_CELL) {
				loaded = 1.0f;
			}
			overhead = ceiling_reach + (std::max(reach, static_cast<float>(h)) - ceiling_reach) * loaded;
		}
		if (p.y - cj * h < overhead && solid(ci, cj - 1)) {
			fy = cj * h;
		} else if ((cj + 1) * h - p.y < reach && solid(ci, cj + 1)) {
			fy = (cj + 1) * h;
		}
		bool hx = !std::isnan(fx);
		bool hy = !std::isnan(fy);
		if (hx) {
			instances.push_back({ Vector2(2.0 * fx - p.x, p.y), Vector2(-v.x, v.y), p_churn[i] });
		}
		if (hy) {
			instances.push_back({ Vector2(p.x, 2.0 * fy - p.y), Vector2(v.x, -v.y), p_churn[i] });
		}
		if (hx && hy) {
			instances.push_back({ Vector2(2.0 * fx - p.x, 2.0 * fy - p.y), -v, p_churn[i] });
		}
	}
	// Against a ceiling the water's top cell is often sparse (one or two particles of a full cell's four, sitting low),
	// and now and then a cell-high layer of air is held there; either way it drew as a dark scalloped strip under every
	// block in a flooded cave. So a cell under a ceiling with water under it is drawn full, as a lattice of splats like
	// a full cell's, mirrored across the faces beside it: when it is water itself, or air in a run that is one cell high
	// and closed at both ends (by blocks or water). A run of air open at an end is real air (a surface with open air
	// beside it), a taller air bell has air under its top row, and a drip under a block has air under it, so none is
	// filled. Water against a side wall is filled the same way, with water beyond it and water or a block over it: never
	// at the water's surface, where filling would draw the water climbing the wall. A mesh of the churned particles
	// alone gets none: the fillers have no churn.
	if (reach > 0.0f && min_churn <= 0.0f) {
		auto air = [&](int i, int j) {
			return _is_cell_in_bounds(i, j) && s[i + f_num_x * j] != 0.0 && cell_types[i + f_num_x * j] == AIR_CELL;
		};
		auto water = [&](int i, int j) {
			return _is_cell_in_bounds(i, j) && cell_types[i + f_num_x * j] == FLUID_CELL;
		};
		auto under_ceiling = [&](int i, int j) { return (air(i, j) || water(i, j)) && solid(i, j - 1) && water(i, j + 1); };
		int per = std::max(particles_per_side, 1);
		auto fill = [&](int c, int j) {
			for (int a = 0; a < per; a++) {
				for (int b = 0; b < per; b++) {
					Vector2 at((c + (a + 0.5) / per) * h, (j + (b + 0.5) / per) * h);
					instances.push_back({ at, Vector2(), 0.0f });
					if (solid(c, j - 1) && at.y - j * h < reach) {
						instances.push_back({ Vector2(at.x, 2.0 * j * h - at.y), Vector2(), 0.0f });
					}
					if (solid(c - 1, j) && at.x - c * h < reach) {
						instances.push_back({ Vector2(2.0 * c * h - at.x, at.y), Vector2(), 0.0f });
					} else if (solid(c + 1, j) && (c + 1) * h - at.x < reach) {
						instances.push_back({ Vector2(2.0 * (c + 1) * h - at.x, at.y), Vector2(), 0.0f });
					}
				}
			}
		};
		for (int j = 1; j < f_num_y - 1; j++) {
			int i = 1;
			while (i < f_num_x - 1) {
				if (!under_ceiling(i, j)) {
					i++;
					continue;
				}
				int start = i;
				bool has_air = false;
				while (i < f_num_x - 1 && under_ceiling(i, j)) {
					has_air = has_air || air(i, j);
					i++;
				}
				bool closed = !air(start - 1, j) && !air(i, j);
				for (int c = start; c < i; c++) {
					if (water(c, j) || (has_air && closed)) {
						fill(c, j);
					}
				}
			}
		}
		for (int j = 1; j < f_num_y - 1; j++) {
			for (int i = 1; i < f_num_x - 1; i++) {
				if (!water(i, j) || under_ceiling(i, j) || !(water(i, j - 1) || solid(i, j - 1))) {
					continue;
				}
				if ((solid(i - 1, j) && water(i + 1, j)) || (solid(i + 1, j) && water(i - 1, j))) {
					fill(i, j);
				}
			}
		}
	}
	int count = static_cast<int>(instances.size());

	// Grow the instance buffer geometrically, so it isn't reallocated every frame. The unused tail is hidden by
	// leaving its transforms zero (a quad shrunk to a point draws nothing), never by lowering the visible instance
	// count: Godot's 2D GL renderer uses the number of instances drawn as the step of its per-item vertex attributes,
	// so every new count is a new vertex layout for the driver to build. A count that follows the particles changes
	// nearly every frame; the layouts piled up by the thousand, and a few minutes in the driver's cache of them
	// overflowed and the AMD driver crashed in its draw thread. All of the buffer is drawn, so the count only
	// changes when the buffer grows.
	int capacity = p_multimesh->get_instance_count();
	if (count > capacity) {
		capacity = std::max(count, capacity + capacity / 2 + 64);
		p_multimesh->set_instance_count(capacity);
	}

	PackedFloat32Array buffer;
	buffer.resize(static_cast<int64_t>(capacity) * stride);
	float *w = buffer.ptrw();
	std::fill(w, w + buffer.size(), 0.0f);
	float scale = static_cast<float>(p_scale);
	for (int n = 0; n < count; n++) {
		const Instance &in = instances[n];
		// 2D layout: basis.x.x, basis.y.x, padding, origin.x, basis.x.y, basis.y.y, padding, origin.y.
		float *t = w + static_cast<int64_t>(n) * stride;
		// Lengthen along the velocity by `stretch` per px/s, up to `max_stretch` times the scale.
		double speed = in.vel.length();
		bool lying = in.lies > 0.75f;
		float along = lying ? scale : scale * static_cast<float>(std::min(1.0 + p_stretch * speed, std::max(p_max_stretch, 1.0)));
		Vector2 dir = speed > 0.001 && !lying ? in.vel / speed : Vector2(1, 0);
		t[0] = dir.x * along;
		t[1] = -dir.y * scale;
		t[3] = in.pos.x;
		t[4] = dir.y * along;
		t[5] = dir.x * scale;
		t[7] = in.pos.y;
		if (colors && p_ground > 0.0) {
			t[8] = in.churn;
			t[9] = in.lies;
			t[10] = lying ? in.floor / static_cast<float>(p_ground) : static_cast<float>(std::min(speed / 600.0, 1.0));
			t[11] = 1.0f;
		} else if (colors) {
			t[8] = t[9] = t[10] = in.churn;
			t[11] = 1.0f;
		}
	}
	p_multimesh->set_buffer(buffer);
	if (p_multimesh->get_visible_instance_count() != -1) {
		p_multimesh->set_visible_instance_count(-1);
	}
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
	q_valid = false;

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
	// Anything that built the query index mid-step saw particles that have since moved.
	q_valid = false;

	frame_effectors.clear();
	volume_sources.clear();
	_update_sensors();
	frame_effector_ids.clear();
	_update_render();
	if (debug_draw) {
		queue_redraw();
	}
}

void FLIPWaterSimulation::_step(double dt) {
	step_dt = dt;
	_apply_effectors(dt);
	_level_thin_water(dt);
	_integrate_particles(dt);
	_push_particles_apart();
	_handle_collisions();
	_update_particle_density();
	_particles_to_grid();
	_solve_incompressibility(dt);
	_apply_wall_drag(dt);
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

// Thin water runs off. A cell under the spray density is not solved: its particles are free, so nothing pushes
// them sideways. That is right for spray in the air, but a thin layer lying on water or on ground just sat where
// it was: the last row of a pool stood above the rim it should have run over, and a film stopped short of a
// ledge's edge. So each row of cells is read for such layers: free water with water or ground under it, carried on
// across solved surface water in the same row and across bare ground (a rim's top, a ledge). Where a run of those
// ends at a drop (air with nothing under it), the layer's particles are pushed toward that end, as their weight
// would push them, up to a walking pace; with a drop at both ends each goes to the nearer. A run shut in at both
// ends (by blocks, by deeper water, by open water at its own level) has nowhere to go and is left alone, so a
// pool's own top row still settles flat. Read from the last step's cells.
void FLIPWaterSimulation::_level_thin_water(double dt) {
	double g = std::abs(gravity.y);
	if (spray_density <= 0.0 || g <= 0.0 || p_positions.empty()) {
		return;
	}
	const double pull = 0.3;
	const double top_speed = 120.0;

	// Particles per cell; thin_push then holds, per cell, which way its free water is to go (0 for nowhere).
	thin_count.assign(f_num_cells, 0);
	for (size_t i = 0; i < p_positions.size(); i++) {
		int xi = std::clamp(static_cast<int>(std::floor(p_positions[i].x * f_inv_spacing)), 0, f_num_x - 1);
		int yi = std::clamp(static_cast<int>(std::floor(p_positions[i].y * f_inv_spacing)), 0, f_num_y - 1);
		thin_count[xi + (f_num_x * yi)]++;
	}
	std::fill(thin_push.begin(), thin_push.end(), 0.0f);

	enum { SHUT, LAYER, PASS, DROP };
	auto kind = [&](int i, int j) {
		int idx = i + (f_num_x * j);
		int under = cell_types[idx + f_num_x];
		if (cell_types[idx] == SOLID_CELL) {
			return SHUT;
		}
		if (cell_types[idx] == FLUID_CELL) {
			// Solved water at the surface carries a layer on; deeper water is a wall to it.
			return cell_types[idx - f_num_x] == AIR_CELL && under != AIR_CELL ? PASS : SHUT;
		}
		if (under == AIR_CELL) {
			return DROP;
		}
		if (thin_count[idx] > 0) {
			// A puddle on the floor stays where it is unless ground run-off is on: it stops the run like a wall.
			return under == SOLID_CELL && !thin_runoff_ground ? SHUT : LAYER;
		}
		return under == SOLID_CELL ? PASS : SHUT;
	};
	bool any = false;
	for (int j = 1; j < f_num_y - 1; j++) {
		int i = 1;
		while (i < f_num_x - 1) {
			int k = kind(i, j);
			if (k != LAYER && k != PASS) {
				i++;
				continue;
			}
			int start = i;
			bool layer = false;
			while (i < f_num_x - 1 && (k = kind(i, j), k == LAYER || k == PASS)) {
				layer = layer || k == LAYER;
				i++;
			}
			bool left = start > 1 && kind(start - 1, j) == DROP;
			bool right = i < f_num_x - 1 && kind(i, j) == DROP;
			if (!layer || (!left && !right)) {
				continue;
			}
			for (int c = start; c < i; c++) {
				if (kind(c, j) != LAYER) {
					continue;
				}
				bool go_left = left && (!right || c - start < i - 1 - c);
				thin_push[c + (f_num_x * j)] = go_left ? -1.0f : 1.0f;
				any = true;
			}
		}
	}
	if (!any) {
		return;
	}
	double gain = pull * g * dt;
	for (size_t i = 0; i < p_positions.size(); i++) {
		int xi = std::clamp(static_cast<int>(std::floor(p_positions[i].x * f_inv_spacing)), 0, f_num_x - 1);
		int yi = std::clamp(static_cast<int>(std::floor(p_positions[i].y * f_inv_spacing)), 0, f_num_y - 1);
		double way = thin_push[xi + (f_num_x * yi)];
		if (way != 0.0 && p_velocities[i].x * way < top_speed) {
			p_velocities[i].x += way * gain;
		}
	}
}

void FLIPWaterSimulation::_integrate_particles(double dt) {
	// Never move more than one cell per step, so particles can't tunnel through walls.
	double limit = cell_size / dt;
	if (max_speed > 0.0) {
		limit = std::min(limit, max_speed);
	}

	// Gravity is already in the velocities: each step ends by adding the next one's (see _grid_to_particles_range).
	for (size_t i = 0; i < p_positions.size(); i++) {
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
		p_churn[kept] = p_churn[i];
		kept++;
	}

	p_positions.resize(kept);
	p_velocities.resize(kept);
	p_churn.resize(kept);
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

	// Churn reads the same neighbour lists, so it runs here, before collisions compact and reorder particles.
	// Cohesion and viscosity share its pass over the neighbours, read from the same snapshot, and land after
	// separation.
	bool cohere = cohesion > 0.0 || viscosity > 0.0;
	if (cohere) {
		p_cohere.resize(num_particles);
		p_smooth.resize(num_particles);
	}
	p_churn.resize(num_particles, 0.0f);
	p_churn_next.resize(num_particles);
	_run_parallel(num_particles, 1024, &FLIPWaterSimulation::_churn_range);
	p_churn.swap(p_churn_next);

	// Jacobi-style: each particle computes its own correction from a read-only snapshot, so this parallelises.
	for (int iter = 0; iter < separation_iterations; iter++) {
		_run_parallel(num_particles, 1024, &FLIPWaterSimulation::_push_particles_range);
		for (int i = 0; i < num_particles; i++) {
			p_positions[i] += p_push[i];
		}
	}

	if (cohere) {
		for (int i = 0; i < num_particles; i++) {
			p_positions[i] += p_cohere[i];
			p_velocities[i] += p_smooth[i];
		}
	}
}

// Churn marks where water traps air: particles slamming into each other at speed, as when a splash lands,
// a jet plunges into a pool or a wave falls onto the water. It is the average speed at which a particle and
// its neighbours (within about a cell) close on each other, mapped from churn_low..churn_high to 0..1. Water
// moving fast but smoothly has neighbours moving with it, so it scores nothing. A particle keeps its peak and
// lets it fade at churn_fade per second, so white water lingers after the impact.
void FLIPWaterSimulation::_churn_range(int begin, int end) {
	double reach = 2.0 * (2.2 * p_radius);
	double reach2 = reach * reach;
	float keep = static_cast<float>(std::exp(-churn_fade * step_dt));
	double span = std::max(churn_high - churn_low, 1.0);
	// Cohesion and viscosity reach two rest spacings; the pull is strongest halfway between one and two.
	bool cohere = cohesion > 0.0 || viscosity > 0.0;
	double rest = 2.0 * p_radius;
	double c_reach = 2.0 * rest;
	double c_reach2 = c_reach * c_reach;
	double c_span = c_reach - rest;

	for (int i = begin; i < end; i++) {
		Vector2 xi = p_positions[i];
		Vector2 vi = p_velocities[i];
		int pxi = static_cast<int>(std::floor(xi.x * p_inv_spacing));
		int pyi = static_cast<int>(std::floor(xi.y * p_inv_spacing));
		double closing = 0.0;
		double weight = 0.0;
		Vector2 pull;
		Vector2 flow;
		double flow_w = 0.0;
		int near = 0;

		for (int xc = std::max(pxi - 2, 0); xc <= std::min(pxi + 2, p_num_x - 1); xc++) {
			for (int yc = std::max(pyi - 2, 0); yc <= std::min(pyi + 2, p_num_y - 1); yc++) {
				int cell_nr = xc + (p_num_x * yc);
				for (int j = first_cell_particle[cell_nr]; j < first_cell_particle[cell_nr + 1]; j++) {
					int id = cell_particle_ids[j];
					if (id == i) {
						continue;
					}
					Vector2 dx = xi - p_positions[id];
					double d2 = dx.length_squared();
					if (d2 > reach2 || d2 == 0.0) {
						continue;
					}
					double d = std::sqrt(d2);
					double w = 1.0 - (d / reach);
					// Positive when the pair is getting closer.
					double approach = -(vi - p_velocities[id]).dot(dx) / d;
					closing += w * std::max(approach, 0.0);
					weight += w;

					if (!cohere || d2 >= c_reach2) {
						continue;
					}
					near++;
					double cw = 1.0 - (d / c_reach);
					flow += (p_velocities[id] - vi) * cw;
					flow_w += cw;
					if (d > rest) {
						double k = (d - rest) / c_span;
						double bump = 4.0 * k * (1.0 - k);
						pull -= dx * (bump * 0.5 * (d - rest) / d);
					}
				}
			}
		}

		// A block beside the particle counts as water at rest: the particle sees, half a spacing inside the block, the
		// neighbours it would have there (one level with it, one a spacing above and one below). Without them a
		// particle at a wall has neighbours on one side only, the pull draws it about half a spacing off the wall, and
		// the surface dips where it meets the wall like a bead of mercury, which reads as water at a small scale.
		// They sit level with the particle, so they pull straight at the wall and never along it; and they only ever
		// cancel a pull away from the wall, never add one toward it, so water does not cling to a wall or creep along
		// it either: a wall is neutral. Floors need none (weight holds water down) and ceilings get none (water must
		// not cling under a block).
		if (cohesion > 0.0) {
			int ci = static_cast<int>(std::floor(xi.x * f_inv_spacing));
			int cj = static_cast<int>(std::floor(xi.y * f_inv_spacing));
			for (int side = -1; side <= 1; side += 2) {
				if (_is_cell_in_bounds(ci + side, cj) && s[(ci + side) + (f_num_x * cj)] != 0.0) {
					continue;
				}
				double face = (side < 0 ? ci : ci + 1) * cell_size;
				double gx = std::abs(xi.x - face) + p_radius;
				double toward = 0.0;
				for (int row = -1; row <= 1; row++) {
					double gy = row * rest;
					double d2 = (gx * gx) + (gy * gy);
					if (d2 >= c_reach2) {
						continue;
					}
					near++;
					double d = std::sqrt(d2);
					if (d > rest) {
						double k = (d - rest) / c_span;
						toward += gx * (4.0 * k * (1.0 - k) * 0.5 * (d - rest) / d);
					}
				}
				double away = -side * pull.x;
				if (away > 0.0) {
					pull.x += side * std::min(toward, away);
				}
			}
		}

		float now = weight > 0.0 ? static_cast<float>(std::clamp(((closing / weight) - churn_low) / span, 0.0, 1.0)) : 0.0f;
		p_churn_next[i] = std::max(p_churn[i] * keep, now);
		if (cohere) {
			p_cohere[i] = near > 0 ? pull * (cohesion / near) : Vector2();
			p_smooth[i] = flow_w > 0.0 ? flow * (viscosity / flow_w) : Vector2();
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
		if (cell_types[cell_nr] == AIR_CELL && (spray_density <= 0.0 || p_densities[cell_nr] >= spray_density * _full_density(cell_nr))) {
			cell_types[cell_nr] = FLUID_CELL;
		}
	}

	// A cell that momentarily holds no particle but is enclosed by fluid is still part of the fluid.
	// Treated as air it would be a zero-pressure hole inside the water, which sucks particles in and
	// keeps the surrounding water churning.
	for (int j = 1; j < f_num_y - 1; j++) {
		for (int i = 1; i < f_num_x - 1; i++) {
			int idx = i + (f_num_x * j);
			if (cell_types[idx] != AIR_CELL) {
				continue;
			}
			const int neighbours[4] = { idx - 1, idx + 1, idx - f_num_x, idx + f_num_x };
			bool enclosed = true;
			bool touches_fluid = false;
			for (int n : neighbours) {
				if (cell_types[n] == AIR_CELL) {
					enclosed = false;
					break;
				}
				touches_fluid = touches_fluid || cell_types[n] == FLUID_CELL;
			}
			if (enclosed && touches_fluid) {
				cell_types[idx] = FLUID_CELL;
			}
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
			const std::vector<double> &wrote = component == 0 ? du : dv;

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
				// A face no particle wrote to holds no flow: reading its zero would drag spray to a halt in mid-air.
				bool open = cell_types[nr[k]] != SOLID_CELL && cell_types[nr[k] - offset] != SOLID_CELL;
				bool solved = cell_types[nr[k]] == FLUID_CELL || cell_types[nr[k] - offset] == FLUID_CELL;
				if (open && (solved || wrote[nr[k]] > 0.0)) {
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
		// The next step's gravity, added now and not as that step starts, so that between steps a particle holds
		// the velocity it is about to move with. Added at the start, water at rest held minus one step of gravity
		// (the solve's answer to the fall it was about to be given), and everything that reads velocities (a body
		// drifting with the current, a splat's stretch) saw still water rising at 16 px/s.
		p_velocities[i] += gravity * step_dt;
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

}

// What a cell's density reads when it is full of water: the rest density, less the share of the density kernel
// that reaches into neighbouring rock (see _compute_drift_source).
double FLIPWaterSimulation::_full_density(int idx) const {
	int x = idx % f_num_x;
	int y = idx / f_num_x;
	if (x == 0 || y == 0 || x == f_num_x - 1 || y == f_num_y - 1) {
		return p_rest_density;
	}
	double a = 0.5 - (p_radius / cell_size);
	double per_rock_side = 1.0 - ((1.0 - a) * (1.0 - a) * 0.5);
	double fx = 1.0;
	double fy = 1.0;
	if (s[idx - 1] == 0.0) {
		fx *= per_rock_side;
	}
	if (s[idx + 1] == 0.0) {
		fx *= per_rock_side;
	}
	if (s[idx - f_num_x] == 0.0) {
		fy *= per_rock_side;
	}
	if (s[idx + f_num_x] == 0.0) {
		fy *= per_rock_side;
	}
	return p_rest_density * fx * fy;
}

// Drift compensation pushes over-compressed fluid apart. Water with a free surface can expand into
// the air, but a sealed pocket (no air cell, no open edge, not touching the unsolved border) cannot
// change volume, so a purely outward push has no solution there: the solver never converges and the
// water churns forever. In sealed pockets the source is made zero-mean, which only redistributes
// density within the pocket.
void FLIPWaterSimulation::_compute_drift_source() {
	std::fill(drift_source.begin(), drift_source.end(), 0.0);
	if (p_rest_density <= 0.0 || drift_compensation <= 0.0) {
		_apply_volume_sources();
		return;
	}

	// Each cell's density is compared with what full water would read there. Next to rock that is less
	// than the rest density: the density kernel reaches half a cell into the rock, and particles keep
	// p_radius off its face, so a share of (1 - a)^2 / 2 of the kernel's weight per rock side stays
	// empty (a = 0.5 - p_radius / h). Comparing wall cells with the plain rest density made full water
	// there read thin, so the pull had to skip them, and water under a block's underside could thin out
	// with nothing pulling it back.
	//
	// Compressed water is pushed apart everywhere. Thinned water is pulled together unless there is air
	// in its 3x3 neighbourhood, since a cell next to the surface reads thin however full it is. Without
	// the pull, pushing alone inflates a body into a sponge of thin pockets; without the push, gravity
	// squeezes it.
	double drift = drift_compensation * cell_size;
	for (int i = 0; i < f_num_cells; i++) {
		if (cell_types[i] != FLUID_CELL) {
			continue;
		}
		int x = i % f_num_x;
		int y = i / f_num_x;
		if (x == 0 || y == 0 || x == f_num_x - 1 || y == f_num_y - 1) {
			double compression = p_densities[i] - p_rest_density;
			if (compression > 0.0) {
				drift_source[i] = drift * compression;
			}
			continue;
		}

		double compression = p_densities[i] - _full_density(i);
		if (compression > 0.0) {
			drift_source[i] = drift * compression;
			continue;
		}
		if (!drift_pull) {
			continue;
		}

		bool near_air = false;
		for (int dy = -1; dy <= 1 && !near_air; dy++) {
			for (int dx = -1; dx <= 1; dx++) {
				if (cell_types[i + dx + (dy * f_num_x)] == AIR_CELL) {
					near_air = true;
					break;
				}
			}
		}
		if (!near_air) {
			drift_source[i] = drift * compression;
		}
	}
	// Before the sealed-pocket pass, so a pocket that cannot change volume evens them out too.
	_apply_volume_sources();

	std::fill(region_ids.begin(), region_ids.end(), -1);
	int region = 0;
	for (int start = 0; start < f_num_cells; start++) {
		if (cell_types[start] != FLUID_CELL || region_ids[start] != -1) {
			continue;
		}

		bool open = false;
		double total = 0.0;
		region_cells.clear();
		region_stack.clear();
		region_stack.push_back(start);
		region_ids[start] = region;
		while (!region_stack.empty()) {
			int idx = region_stack.back();
			region_stack.pop_back();
			region_cells.push_back(idx);
			total += drift_source[idx];

			int i = idx % f_num_x;
			int j = idx / f_num_x;
			// Border cells are never solved, so they act as a free (zero pressure) boundary.
			if (i == 0 || j == 0 || i == f_num_x - 1 || j == f_num_y - 1) {
				open = true;
				continue;
			}
			const int neighbours[4] = { idx - 1, idx + 1, idx - f_num_x, idx + f_num_x };
			for (int n : neighbours) {
				if (cell_types[n] == AIR_CELL) {
					open = true;
				} else if (cell_types[n] == FLUID_CELL && region_ids[n] == -1) {
					region_ids[n] = region;
					region_stack.push_back(n);
				}
			}
		}

		if (!open && !region_cells.empty()) {
			double mean = total / region_cells.size();
			for (int idx : region_cells) {
				drift_source[idx] -= mean;
			}
		}
		region++;
	}
}

// A swirl flings particles out of its core, and the emptied cells would count as open air at zero pressure: a hole
// opens in the water in a frame and slams shut when it collapses. Air enclosed by water (not reaching the border, not
// held under a solid) and no bigger than void_fill_cells is solved as water instead. It reads empty, so drift
// compensation pulls particles in and the core thins and refills gradually. Air bells under blocks, and big pockets
// like a breaking wave's tube, stay air.
void FLIPWaterSimulation::_fill_voids() {
	if (void_fill_cells <= 0) {
		return;
	}
	std::fill(region_ids.begin(), region_ids.end(), -1);
	for (int start = 0; start < f_num_cells; start++) {
		if (cell_types[start] != AIR_CELL || region_ids[start] != -1) {
			continue;
		}
		bool sealed = true;
		region_cells.clear();
		region_stack.clear();
		region_stack.push_back(start);
		region_ids[start] = 0;
		while (!region_stack.empty()) {
			int idx = region_stack.back();
			region_stack.pop_back();
			region_cells.push_back(idx);
			int i = idx % f_num_x;
			int j = idx / f_num_x;
			if (i == 0 || j == 0 || i == f_num_x - 1 || j == f_num_y - 1) {
				sealed = false;
				continue;
			}
			const int neighbours[4] = { idx - 1, idx + 1, idx - f_num_x, idx + f_num_x };
			// Air under a ceiling is a held air bell; a block beside or below the pocket doesn't hold it.
			if (cell_types[idx - f_num_x] == SOLID_CELL) {
				sealed = false;
			}
			for (int nb : neighbours) {
				if (cell_types[nb] == AIR_CELL && region_ids[nb] == -1) {
					region_ids[nb] = 0;
					region_stack.push_back(nb);
				}
			}
		}
		if (sealed && static_cast<int>(region_cells.size()) <= void_fill_cells) {
			for (int idx : region_cells) {
				cell_types[idx] = FLUID_CELL;
			}
		}
	}
}

void FLIPWaterSimulation::_solve_incompressibility(double dt) {
	std::fill(pressure.begin(), pressure.end(), 0.0);
	prev_u = u;
	prev_v = v;
	solve_cp = DENSITY * cell_size / dt;
	_fill_voids();
	_compute_drift_source();

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

// Floors are otherwise frictionless: a jet hitting a ledge would turn all its speed sideways and fly off as a sheet.
// Flow skimming over a solid below it (u faces in a cell on a floor) is braked like turbulent friction, by a rate that
// grows with its speed: wall_drag per second at 500 px/s, a fifth of that at 100 px/s. Fast sheets lose their speed in
// a block or two while slow water still runs freely, and walls and ceilings stay frictionless so nothing clings.
void FLIPWaterSimulation::_apply_wall_drag(double dt) {
	if (wall_drag <= 0.0) {
		return;
	}
	double c = wall_drag * dt / 500.0;
	auto solid = [&](int idx) { return cell_types[idx] == SOLID_CELL; };
	for (int j = 1; j < f_num_y - 1; j++) {
		for (int i = 1; i < f_num_x - 1; i++) {
			int idx = i + (f_num_x * j);
			if (solid(idx) || solid(idx - 1) || (cell_types[idx] != FLUID_CELL && cell_types[idx - 1] != FLUID_CELL)) {
				continue;
			}
			int down = idx + f_num_x;
			int floor = solid(down) ? down : solid(down - 1) ? down - 1 : -1;
			if (floor < 0) {
				continue;
			}
			// Implicit quadratic drag, so it never overshoots past the floor's own speed.
			double target = solid_velocity[floor].x;
			double rel = u[idx] - target;
			u[idx] = target + rel / (1.0 + c * std::abs(rel));
		}
	}
}

void FLIPWaterSimulation::_solve_rows_range(int begin, int end) {
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

			double div = u[right] - u[center] + v[bottom] - v[center] - drift_source[center];

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
