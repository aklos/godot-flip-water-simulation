#ifndef FLIP_EFFECTORS_H
#define FLIP_EFFECTORS_H

#include <godot_cpp/classes/node2d.hpp>

namespace godot {

class FLIPWaterSimulation;

// Base class for nodes that act on a FLIPWaterSimulation. The simulation is found through
// `simulation`, or the nearest FLIPWaterSimulation ancestor when that path is empty.
class FLIPEffector : public Node2D {
	GDCLASS(FLIPEffector, Node2D)

	NodePath simulation;
	bool enabled = true;

protected:
	static void _bind_methods();
	void _notification(int p_what);
	virtual void _draw_gizmo() {}

public:
	static constexpr const char *GROUP = "_flip_effectors";

	void set_simulation(const NodePath &p_path);
	NodePath get_simulation() const;
	void set_enabled(bool p_enabled);
	bool is_enabled() const;

	FLIPWaterSimulation *get_simulation_node() const;
};

// Spawns particles along a line of `width`, moving along the node's +X axis.
class FLIPEmitter : public FLIPEffector {
	GDCLASS(FLIPEmitter, FLIPEffector)

	double rate = 100.0;
	double speed = 300.0;
	double spread = 0.2;
	double width = 20.0;
	double accumulator = 0.0;
	int pending_burst = 0;

protected:
	static void _bind_methods();
	void _draw_gizmo() override;

public:
	void set_rate(double p_rate);
	double get_rate() const;
	void set_speed(double p_speed);
	double get_speed() const;
	void set_spread(double p_spread);
	double get_spread() const;
	void set_width(double p_width);
	double get_width() const;
	void emit_burst(int p_count);

	int take_spawn_count(double p_dt);
};

// Removes particles within `radius`.
class FLIPDrain : public FLIPEffector {
	GDCLASS(FLIPDrain, FLIPEffector)

	double radius = 30.0;
	int drained_this_frame = 0;

protected:
	static void _bind_methods();
	void _draw_gizmo() override;

public:
	void set_radius(double p_radius);
	double get_radius() const;

	void add_drained(int p_count);
	void flush_drained();
};

// Accelerates particles within `radius`, either away from the node (RADIAL) or along its +X axis
// (DIRECTIONAL). `spread` limits the effect to a cone around +X.
class FLIPForce : public FLIPEffector {
	GDCLASS(FLIPForce, FLIPEffector)

public:
	enum ForceMode {
		MODE_RADIAL = 0,
		MODE_DIRECTIONAL = 1,
	};

private:
	ForceMode mode = MODE_RADIAL;
	double strength = 2000.0;
	double radius = 100.0;
	double spread = Math_TAU;
	bool falloff = true;

protected:
	static void _bind_methods();
	void _draw_gizmo() override;

public:
	void set_mode(ForceMode p_mode);
	ForceMode get_mode() const;
	void set_strength(double p_strength);
	double get_strength() const;
	void set_radius(double p_radius);
	double get_radius() const;
	void set_spread(double p_spread);
	double get_spread() const;
	void set_falloff(bool p_falloff);
	bool has_falloff() const;
};

// Counts particles inside a rectangle of `size` centred on the node and signals when water arrives or leaves.
class FLIPSensor : public FLIPEffector {
	GDCLASS(FLIPSensor, FLIPEffector)

	Vector2 size = Vector2(64, 64);
	int threshold = 1;
	int particle_count = 0;

protected:
	static void _bind_methods();
	void _draw_gizmo() override;

public:
	void set_size(const Vector2 &p_size);
	Vector2 get_size() const;
	void set_threshold(int p_threshold);
	int get_threshold() const;
	int get_particle_count() const;
	bool is_wet() const;

	void set_particle_count(int p_count);
};

} // namespace godot

VARIANT_ENUM_CAST(godot::FLIPForce::ForceMode);

#endif // FLIP_EFFECTORS_H
