#include "flip_effectors.h"
#include "flip_sim.h"

#include <godot_cpp/classes/engine.hpp>
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/math.hpp>

#include <algorithm>
#include <cmath>

using namespace godot;

static const Color GIZMO_COLOR = Color(0.35, 0.65, 1.0, 0.8);

// FLIPEffector

void FLIPEffector::_bind_methods() {
	ClassDB::bind_method(D_METHOD("set_simulation", "path"), &FLIPEffector::set_simulation);
	ClassDB::bind_method(D_METHOD("get_simulation"), &FLIPEffector::get_simulation);
	ClassDB::bind_method(D_METHOD("set_enabled", "enabled"), &FLIPEffector::set_enabled);
	ClassDB::bind_method(D_METHOD("is_enabled"), &FLIPEffector::is_enabled);

	ADD_PROPERTY(PropertyInfo(Variant::NODE_PATH, "simulation", PROPERTY_HINT_NODE_PATH_VALID_TYPES, "FLIPWaterSimulation"), "set_simulation", "get_simulation");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "enabled"), "set_enabled", "is_enabled");
}

void FLIPEffector::_notification(int p_what) {
	switch (p_what) {
		case NOTIFICATION_ENTER_TREE:
			add_to_group(GROUP);
			break;
		case NOTIFICATION_DRAW:
			if (Engine::get_singleton()->is_editor_hint()) {
				_draw_gizmo();
			}
			break;
	}
}

void FLIPEffector::set_simulation(const NodePath &p_path) {
	simulation = p_path;
}

NodePath FLIPEffector::get_simulation() const {
	return simulation;
}

void FLIPEffector::set_enabled(bool p_enabled) {
	enabled = p_enabled;
}

bool FLIPEffector::is_enabled() const {
	return enabled;
}

FLIPWaterSimulation *FLIPEffector::get_simulation_node() const {
	if (!is_inside_tree()) {
		return nullptr;
	}
	if (!simulation.is_empty()) {
		return Object::cast_to<FLIPWaterSimulation>(get_node_or_null(simulation));
	}
	for (Node *n = get_parent(); n; n = n->get_parent()) {
		if (FLIPWaterSimulation *sim = Object::cast_to<FLIPWaterSimulation>(n)) {
			return sim;
		}
	}
	return nullptr;
}

// FLIPEmitter

void FLIPEmitter::_bind_methods() {
	ClassDB::bind_method(D_METHOD("set_rate", "rate"), &FLIPEmitter::set_rate);
	ClassDB::bind_method(D_METHOD("get_rate"), &FLIPEmitter::get_rate);
	ClassDB::bind_method(D_METHOD("set_speed", "speed"), &FLIPEmitter::set_speed);
	ClassDB::bind_method(D_METHOD("get_speed"), &FLIPEmitter::get_speed);
	ClassDB::bind_method(D_METHOD("set_spread", "spread"), &FLIPEmitter::set_spread);
	ClassDB::bind_method(D_METHOD("get_spread"), &FLIPEmitter::get_spread);
	ClassDB::bind_method(D_METHOD("set_width", "width"), &FLIPEmitter::set_width);
	ClassDB::bind_method(D_METHOD("get_width"), &FLIPEmitter::get_width);
	ClassDB::bind_method(D_METHOD("emit_burst", "count"), &FLIPEmitter::emit_burst);

	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "rate", PROPERTY_HINT_RANGE, "0,5000,1,or_greater,suffix:/s"), "set_rate", "get_rate");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "speed", PROPERTY_HINT_RANGE, "0,2000,1,or_greater,suffix:px/s"), "set_speed", "get_speed");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "spread", PROPERTY_HINT_RANGE, "0,360,0.1,radians_as_degrees"), "set_spread", "get_spread");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "width", PROPERTY_HINT_RANGE, "0,500,0.1,or_greater,suffix:px"), "set_width", "get_width");
}

void FLIPEmitter::_draw_gizmo() {
	draw_line(Vector2(0, -width * 0.5), Vector2(0, width * 0.5), GIZMO_COLOR, 2.0);
	double len = std::clamp(speed * 0.1, 16.0, 128.0);
	draw_line(Vector2(), Vector2(len, 0), GIZMO_COLOR, 2.0);
	draw_line(Vector2(len, 0), Vector2(len - 6, -4), GIZMO_COLOR, 2.0);
	draw_line(Vector2(len, 0), Vector2(len - 6, 4), GIZMO_COLOR, 2.0);
}

void FLIPEmitter::set_rate(double p_rate) {
	rate = std::max(p_rate, 0.0);
}

double FLIPEmitter::get_rate() const {
	return rate;
}

void FLIPEmitter::set_speed(double p_speed) {
	speed = p_speed;
	queue_redraw();
}

double FLIPEmitter::get_speed() const {
	return speed;
}

void FLIPEmitter::set_spread(double p_spread) {
	spread = std::clamp(p_spread, 0.0, Math_TAU);
}

double FLIPEmitter::get_spread() const {
	return spread;
}

void FLIPEmitter::set_width(double p_width) {
	width = std::max(p_width, 0.0);
	queue_redraw();
}

double FLIPEmitter::get_width() const {
	return width;
}

void FLIPEmitter::emit_burst(int p_count) {
	pending_burst += std::max(p_count, 0);
}

int FLIPEmitter::take_spawn_count(double p_dt) {
	accumulator += rate * p_dt;
	int count = static_cast<int>(accumulator);
	accumulator -= count;
	count += pending_burst;
	pending_burst = 0;
	return count;
}

// FLIPDrain

void FLIPDrain::_bind_methods() {
	ClassDB::bind_method(D_METHOD("set_radius", "radius"), &FLIPDrain::set_radius);
	ClassDB::bind_method(D_METHOD("get_radius"), &FLIPDrain::get_radius);

	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "radius", PROPERTY_HINT_RANGE, "0,500,0.1,or_greater,suffix:px"), "set_radius", "get_radius");

	ADD_SIGNAL(MethodInfo("water_drained", PropertyInfo(Variant::INT, "count")));
}

void FLIPDrain::_draw_gizmo() {
	draw_arc(Vector2(), radius, 0, Math_TAU, 48, GIZMO_COLOR, 2.0);
}

void FLIPDrain::set_radius(double p_radius) {
	radius = std::max(p_radius, 0.0);
	queue_redraw();
}

double FLIPDrain::get_radius() const {
	return radius;
}

void FLIPDrain::add_drained(int p_count) {
	drained_this_frame += p_count;
}

void FLIPDrain::flush_drained() {
	if (drained_this_frame > 0) {
		int count = drained_this_frame;
		drained_this_frame = 0;
		emit_signal("water_drained", count);
	}
}

// FLIPForce

void FLIPForce::_bind_methods() {
	ClassDB::bind_method(D_METHOD("set_mode", "mode"), &FLIPForce::set_mode);
	ClassDB::bind_method(D_METHOD("get_mode"), &FLIPForce::get_mode);
	ClassDB::bind_method(D_METHOD("set_strength", "strength"), &FLIPForce::set_strength);
	ClassDB::bind_method(D_METHOD("get_strength"), &FLIPForce::get_strength);
	ClassDB::bind_method(D_METHOD("set_radius", "radius"), &FLIPForce::set_radius);
	ClassDB::bind_method(D_METHOD("get_radius"), &FLIPForce::get_radius);
	ClassDB::bind_method(D_METHOD("set_spread", "spread"), &FLIPForce::set_spread);
	ClassDB::bind_method(D_METHOD("get_spread"), &FLIPForce::get_spread);
	ClassDB::bind_method(D_METHOD("set_falloff", "falloff"), &FLIPForce::set_falloff);
	ClassDB::bind_method(D_METHOD("has_falloff"), &FLIPForce::has_falloff);

	ADD_PROPERTY(PropertyInfo(Variant::INT, "mode", PROPERTY_HINT_ENUM, "Radial,Directional"), "set_mode", "get_mode");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "strength", PROPERTY_HINT_RANGE, "-10000,10000,1,or_greater,or_less,suffix:px/s²"), "set_strength", "get_strength");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "radius", PROPERTY_HINT_RANGE, "0,2000,0.1,or_greater,suffix:px"), "set_radius", "get_radius");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "spread", PROPERTY_HINT_RANGE, "0,360,0.1,radians_as_degrees"), "set_spread", "get_spread");
	ADD_PROPERTY(PropertyInfo(Variant::BOOL, "falloff"), "set_falloff", "has_falloff");

	BIND_ENUM_CONSTANT(MODE_RADIAL);
	BIND_ENUM_CONSTANT(MODE_DIRECTIONAL);
}

void FLIPForce::_draw_gizmo() {
	if (spread >= Math_TAU - 0.001) {
		draw_arc(Vector2(), radius, 0, Math_TAU, 48, GIZMO_COLOR, 2.0);
	} else {
		double half = spread * 0.5;
		draw_arc(Vector2(), radius, -half, half, 24, GIZMO_COLOR, 2.0);
		draw_line(Vector2(), Vector2(radius, 0).rotated(-half), GIZMO_COLOR, 2.0);
		draw_line(Vector2(), Vector2(radius, 0).rotated(half), GIZMO_COLOR, 2.0);
	}
	if (mode == MODE_DIRECTIONAL) {
		double len = std::min(radius, 48.0);
		draw_line(Vector2(), Vector2(len, 0), GIZMO_COLOR, 2.0);
		draw_line(Vector2(len, 0), Vector2(len - 6, -4), GIZMO_COLOR, 2.0);
		draw_line(Vector2(len, 0), Vector2(len - 6, 4), GIZMO_COLOR, 2.0);
	}
}

void FLIPForce::set_mode(ForceMode p_mode) {
	mode = p_mode;
	queue_redraw();
}

FLIPForce::ForceMode FLIPForce::get_mode() const {
	return mode;
}

void FLIPForce::set_strength(double p_strength) {
	strength = p_strength;
}

double FLIPForce::get_strength() const {
	return strength;
}

void FLIPForce::set_radius(double p_radius) {
	radius = std::max(p_radius, 0.0);
	queue_redraw();
}

double FLIPForce::get_radius() const {
	return radius;
}

void FLIPForce::set_spread(double p_spread) {
	spread = std::clamp(p_spread, 0.0, Math_TAU);
	queue_redraw();
}

double FLIPForce::get_spread() const {
	return spread;
}

void FLIPForce::set_falloff(bool p_falloff) {
	falloff = p_falloff;
}

bool FLIPForce::has_falloff() const {
	return falloff;
}

// FLIPSensor

void FLIPSensor::_bind_methods() {
	ClassDB::bind_method(D_METHOD("set_size", "size"), &FLIPSensor::set_size);
	ClassDB::bind_method(D_METHOD("get_size"), &FLIPSensor::get_size);
	ClassDB::bind_method(D_METHOD("set_threshold", "threshold"), &FLIPSensor::set_threshold);
	ClassDB::bind_method(D_METHOD("get_threshold"), &FLIPSensor::get_threshold);
	ClassDB::bind_method(D_METHOD("get_particle_count"), &FLIPSensor::get_particle_count);
	ClassDB::bind_method(D_METHOD("is_wet"), &FLIPSensor::is_wet);

	ADD_PROPERTY(PropertyInfo(Variant::VECTOR2, "size", PROPERTY_HINT_NONE, "suffix:px"), "set_size", "get_size");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "threshold", PROPERTY_HINT_RANGE, "1,1000,1,or_greater"), "set_threshold", "get_threshold");

	ADD_SIGNAL(MethodInfo("water_entered"));
	ADD_SIGNAL(MethodInfo("water_exited"));
}

void FLIPSensor::_draw_gizmo() {
	draw_rect(Rect2(-size * 0.5, size), GIZMO_COLOR, false, 2.0);
}

void FLIPSensor::set_size(const Vector2 &p_size) {
	size = p_size.abs();
	queue_redraw();
}

Vector2 FLIPSensor::get_size() const {
	return size;
}

void FLIPSensor::set_threshold(int p_threshold) {
	threshold = std::max(p_threshold, 1);
}

int FLIPSensor::get_threshold() const {
	return threshold;
}

int FLIPSensor::get_particle_count() const {
	return particle_count;
}

bool FLIPSensor::is_wet() const {
	return particle_count >= threshold;
}

void FLIPSensor::set_particle_count(int p_count) {
	bool was_wet = is_wet();
	particle_count = p_count;
	bool wet = is_wet();
	if (wet && !was_wet) {
		emit_signal("water_entered");
	} else if (!wet && was_wet) {
		emit_signal("water_exited");
	}
}
