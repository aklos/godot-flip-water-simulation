#ifndef GDEXAMPLE_H
#define GDEXAMPLE_H

#include <godot_cpp/classes/node2d.hpp>

namespace godot {

class FLIPWaterSimulation : public Node2D {
	GDCLASS(FLIPWaterSimulation, Node2D)

private:
    double cell_size;
    double density;
    int f_num_x;
    int f_num_y;
    int f_num_cells;
    double f_inv_spacing;
    double p_radius;
    double p_rest_density;
    int AIR_CELL;
    int SOLID_CELL;
    int FLUID_CELL;
    double p_inv_spacing;
    int p_num_x;
    int p_num_y;
    std::vector<int> cell_types;
    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> du;
    std::vector<double> dv;
    std::vector<double> prev_u;
    std::vector<double> prev_v;
    std::vector<double> pressure;
    std::vector<double> s;
    std::vector<double> p_densities;
    std::vector<Vector2> p_positions;
    std::vector<Vector2> p_velocities;
    std::vector<int> num_cell_particles;
    std::vector<int> first_cell_particle;
    std::vector<int> cell_particle_ids;
    std::vector<Vector2> jet_positions;
    std::vector<Vector2> jet_directions;
    int _clamp_int(int v, int min, int max);
    double _clamp_double(double v, double min, double max);
    void _integrate_particles(double delta);
    void _handle_collisions();
    void _push_particles_apart();
    void _transfer_velocities(bool to_grid);
    void _update_particle_density();
    void _solve_incompressibility(double delta, int iters);

protected:
	static void _bind_methods();

public:
	FLIPWaterSimulation();
	~FLIPWaterSimulation();

    void init(const double width, const double height, const double scale);
    double get_cell_size();
    double get_particle_radius();
	void update(double delta);
    void shift_down();
    void shift_up();
    void set_row(const int row, const Array r_cell_types);
    void set_cell_types(const Array r_cell_types);
    void add_wall(Vector2 pos);
    void remove_wall(Vector2 pos);
    void push_particle(int index, Vector2 vel);
    void add_particle(Vector2 pos);
    void remove_particle(int index);
    void add_jet(Vector2 pos, Vector2 dir);
    Array get_particle_data();
    Array get_cell_types();
    Array get_pressure();
};

}

#endif