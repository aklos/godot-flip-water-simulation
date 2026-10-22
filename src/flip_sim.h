#ifndef FLIP_SIM_H
#define FLIP_SIM_H

#include <godot_cpp/classes/node2d.hpp>
#include <godot_cpp/classes/physics_shape_query_parameters2d.hpp>
#include <godot_cpp/classes/quad_mesh.hpp>
#include <godot_cpp/classes/shader_material.hpp>
#include <godot_cpp/variant/packed_float32_array.hpp>
#include <godot_cpp/variant/packed_int32_array.hpp>
#include <godot_cpp/variant/packed_vector2_array.hpp>
#include <godot_cpp/variant/typed_array.hpp>

#include <cstdint>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace godot {

class FLIPEffector;
class FLIPEmitter;
class FLIPDrain;
class FLIPForce;
class Sprite2D;
class SubViewport;

class FLIPWaterSimulation : public Node2D {
	GDCLASS(FLIPWaterSimulation, Node2D)

public:
	enum CellType {
		AIR_CELL = 0,
		SOLID_CELL = 1,
		FLUID_CELL = 2,
	};

	enum EdgeMode {
		EDGE_WALL = 0,
		EDGE_OPEN = 1,
	};

	enum RenderMode {
		RENDER_NONE = 0,
		RENDER_PARTICLES = 1,
		RENDER_SURFACE = 2,
	};

private:
	// Grid.
	Vector2i grid_size;
	double cell_size;
	int f_num_x;
	int f_num_y;
	int f_num_cells;
	double f_inv_spacing;
	double p_radius;
	double p_rest_density;
	double p_inv_spacing;
	int p_num_x;
	int p_num_y;

	// Settings.
	Vector2 gravity;
	double flip_ratio;
	double max_speed;
	int solver_iterations;
	int substeps;
	int max_particles;
	EdgeMode edge_modes[4];
	bool auto_update;
	bool use_threads;
	NodePath follow_target;
	uint32_t collision_mask;
	bool debug_draw;

	// Rendering settings.
	RenderMode render_mode;
	Color water_color;
	Color foam_color;
	double foam_speed;
	double particle_scale;
	double surface_resolution;

	// Cell state. A cell is solid if any of the solid layers (or a WALL edge) marks it.
	std::vector<int> cell_types;
	std::vector<uint8_t> user_solid;
	std::vector<uint8_t> static_solid;
	std::vector<uint8_t> dynamic_solid;
	std::vector<Vector2> solid_velocity;
	std::vector<double> u;
	std::vector<double> v;
	std::vector<double> du;
	std::vector<double> dv;
	std::vector<double> prev_u;
	std::vector<double> prev_v;
	std::vector<double> pressure;
	std::vector<double> s;
	std::vector<double> p_densities;

	// Particle state.
	std::vector<Vector2> p_positions;
	std::vector<Vector2> p_velocities;
	std::vector<Vector2> p_push;
	std::vector<int> num_cell_particles;
	std::vector<int> first_cell_particle;
	std::vector<int> cell_particle_ids;
	std::minstd_rand rng;

	// Physics colliders.
	RID query_shape;
	Ref<PhysicsShapeQueryParameters2D> query_params;
	std::unordered_set<uint64_t> static_bodies;
	TypedArray<RID> static_body_rids;
	std::unordered_map<uint64_t, Transform2D> body_prev_transforms;
	std::vector<int> dynamic_cells;
	bool static_colliders_dirty;

	// Multithreading: one job at a time, split into chunks across the WorkerThreadPool.
	void (FLIPWaterSimulation::*job_fn)(int, int);
	int job_count;
	int job_chunks;
	int solve_color;
	double solve_cp;

	// Rendering.
	RID render_item;
	RID render_multimesh;
	Ref<QuadMesh> render_quad;
	Ref<ShaderMaterial> particle_material;
	Ref<ShaderMaterial> splat_material;
	Ref<ShaderMaterial> surface_material;
	SubViewport *surface_viewport;
	Node2D *surface_root;
	Sprite2D *surface_sprite;
	PackedFloat32Array render_buffer;
	int render_capacity;

	// Effectors gathered for the current update(). IDs are re-resolved before emitting signals,
	// since signal handlers may free nodes.
	std::vector<FLIPEffector *> frame_effectors;
	std::vector<uint64_t> frame_effector_ids;

	void _rebuild_grid();
	bool _is_cell_in_bounds(int i, int j) const;
	bool _is_edge_wall(int i, int j) const;
	void _refresh_cell(int i, int j);
	void _update_solids();
	void _spawn_cell_particles(double x, double y);
	bool _add_particle_internal(const Vector2 &pos, const Vector2 &vel);
	void _remove_particles_outside_domain();
	void _follow();

	void _refresh_static_colliders();
	void _update_dynamic_colliders(double delta);
	void _rasterize_shape(const RID &body, int shape_idx, bool dynamic, const Transform2D &sim_xform, const Transform2D &body_xform, double delta);
	void _gather_effectors();
	void _apply_effectors(double dt);
	void _apply_emitter(FLIPEmitter *emitter, const Transform2D &xf, double dt);
	void _apply_force(FLIPForce *force, const Transform2D &xf, double dt);
	int _apply_drain(FLIPDrain *drain, const Transform2D &xf);
	void _update_sensors();

	void _step(double dt);
	void _integrate_particles(double dt);
	void _push_particles_apart();
	void _push_particles_range(int begin, int end);
	void _handle_collisions();
	void _particles_to_grid();
	void _grid_to_particles_range(int begin, int end);
	void _update_particle_density();
	void _solve_incompressibility(double dt);
	void _solve_rows_range(int begin, int end);

	void _run_parallel(int count, int min_per_chunk, void (FLIPWaterSimulation::*fn)(int, int));
	void _run_job_chunk(uint32_t chunk);

	void _ensure_render_resources();
	void _free_render_resources();
	void _update_render();

protected:
	static void _bind_methods();
	void _notification(int p_what);

public:
	FLIPWaterSimulation();
	~FLIPWaterSimulation();

	void set_grid_size(const Vector2i &p_size);
	Vector2i get_grid_size() const;
	void set_cell_size(double p_size);
	double get_cell_size() const;
	void set_gravity(const Vector2 &p_gravity);
	Vector2 get_gravity() const;
	void set_flip_ratio(double p_flip_ratio);
	double get_flip_ratio() const;
	void set_max_speed(double p_max_speed);
	double get_max_speed() const;
	void set_solver_iterations(int p_iterations);
	int get_solver_iterations() const;
	void set_substeps(int p_substeps);
	int get_substeps() const;
	void set_max_particles(int p_max);
	int get_max_particles() const;
	void set_edge_mode(int p_side, EdgeMode p_mode);
	EdgeMode get_edge_mode(int p_side) const;
	void set_auto_update(bool p_enabled);
	bool is_auto_update() const;
	void set_use_threads(bool p_enabled);
	bool is_using_threads() const;
	void set_follow_target(const NodePath &p_path);
	NodePath get_follow_target() const;
	void set_collision_mask(uint32_t p_mask);
	uint32_t get_collision_mask() const;
	void set_debug_draw(bool p_enabled);
	bool is_debug_draw() const;

	void set_render_mode(RenderMode p_mode);
	RenderMode get_render_mode() const;
	void set_water_color(const Color &p_color);
	Color get_water_color() const;
	void set_foam_color(const Color &p_color);
	Color get_foam_color() const;
	void set_foam_speed(double p_speed);
	double get_foam_speed() const;
	void set_particle_scale(double p_scale);
	double get_particle_scale() const;
	void set_surface_resolution(double p_resolution);
	double get_surface_resolution() const;

	double get_particle_radius() const;
	Rect2 get_domain_rect() const;
	Vector2i local_to_cell(const Vector2 &p_pos) const;
	Vector2 cell_to_local(const Vector2i &p_cell) const;

	void set_cell(const Vector2i &p_cell, CellType p_type);
	CellType get_cell(const Vector2i &p_cell) const;
	void fill_cells(const Rect2i &p_rect, CellType p_type);
	void clear();
	void shift(const Vector2i &p_cells);
	void refresh_colliders();
	void update(double delta);

	int get_particle_count() const;
	PackedVector2Array get_particle_positions() const;
	PackedVector2Array get_particle_velocities() const;
	PackedInt32Array get_cell_types() const;
	PackedFloat32Array get_pressure() const;
	bool add_particle(const Vector2 &p_pos, const Vector2 &p_velocity = Vector2());
	void apply_impulse(const Vector2 &p_center, double p_radius, const Vector2 &p_velocity);
	int remove_particles_in_circle(const Vector2 &p_center, double p_radius);
	int count_particles_in_circle(const Vector2 &p_center, double p_radius) const;
	int count_particles_in_rect(const Rect2 &p_rect) const;
};

} // namespace godot

VARIANT_ENUM_CAST(godot::FLIPWaterSimulation::CellType);
VARIANT_ENUM_CAST(godot::FLIPWaterSimulation::EdgeMode);
VARIANT_ENUM_CAST(godot::FLIPWaterSimulation::RenderMode);

#endif // FLIP_SIM_H
