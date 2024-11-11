#include "flip_sim.h"
#include <godot_cpp/core/class_db.hpp>
#include <algorithm>
// #include <execution>

using namespace godot;

void FLIPWaterSimulation::_bind_methods() {
	ClassDB::bind_method(D_METHOD("init", "width", "height", "scale"), &FLIPWaterSimulation::init);
	ClassDB::bind_method(D_METHOD("get_cell_size"), &FLIPWaterSimulation::get_cell_size);
	ClassDB::bind_method(D_METHOD("get_particle_radius"), &FLIPWaterSimulation::get_particle_radius);
	ClassDB::bind_method(D_METHOD("get_particle_data"), &FLIPWaterSimulation::get_particle_data);
	ClassDB::bind_method(D_METHOD("get_cell_types"), &FLIPWaterSimulation::get_cell_types);
	ClassDB::bind_method(D_METHOD("get_pressure"), &FLIPWaterSimulation::get_pressure);
	ClassDB::bind_method(D_METHOD("update", "delta"), &FLIPWaterSimulation::update);
	ClassDB::bind_method(D_METHOD("set_row", "row", "r_cell_types"), &FLIPWaterSimulation::set_row);
	ClassDB::bind_method(D_METHOD("set_cell_types", "r_cell_types"), &FLIPWaterSimulation::set_cell_types);
	ClassDB::bind_method(D_METHOD("add_wall", "pos"), &FLIPWaterSimulation::add_wall);
	ClassDB::bind_method(D_METHOD("remove_wall", "pos"), &FLIPWaterSimulation::remove_wall);
	ClassDB::bind_method(D_METHOD("push_particle", "index", "vel"), &FLIPWaterSimulation::push_particle);
	ClassDB::bind_method(D_METHOD("add_particle", "pos"), &FLIPWaterSimulation::add_particle);
	ClassDB::bind_method(D_METHOD("remove_particle", "index"), &FLIPWaterSimulation::remove_particle);
	ClassDB::bind_method(D_METHOD("add_jet", "pos", "dir"), &FLIPWaterSimulation::add_jet);
	ClassDB::bind_method(D_METHOD("shift_down"), &FLIPWaterSimulation::shift_down);
	ClassDB::bind_method(D_METHOD("shift_up"), &FLIPWaterSimulation::shift_up);
}

FLIPWaterSimulation::FLIPWaterSimulation() {
	cell_size = 0.0;
	density = 1000.0;
	f_num_x = 0;
	f_num_y = 0;
	f_num_cells = 0;
	p_rest_density = 0.0;
	p_radius = 0.0;
	f_inv_spacing = 0.0;
	AIR_CELL = 0;
	SOLID_CELL = 1;
	FLUID_CELL = 2;
	p_num_x = 0;
	p_num_y = 0;
	p_inv_spacing = 0.0;
}

FLIPWaterSimulation::~FLIPWaterSimulation() {
	// Add your cleanup here.
}

void FLIPWaterSimulation::init(
	const double width,
	const double height,
	const double scale
) {
	cell_size = 120.0 * scale / 2.0;
	density = 1000.0;
	f_num_x = static_cast<int>(std::floor(width / cell_size));
	f_num_y = static_cast<int>(std::floor(height / cell_size));
	f_num_cells = f_num_x * f_num_y;
	p_rest_density = 0.0;
	p_radius = 0.3 * cell_size;
	f_inv_spacing = 1.0 / cell_size;
	u.resize(f_num_cells, 0.0);
	v.resize(f_num_cells, 0.0);
	du.resize(f_num_cells, 0.0);
	dv.resize(f_num_cells, 0.0);
	prev_u.resize(f_num_cells, 0.0);
	prev_v.resize(f_num_cells, 0.0);
	pressure.resize(f_num_cells, 0.0);
	p_densities.resize(f_num_cells, 0.0);
	cell_types.resize(f_num_cells, 0);
	s.resize(f_num_cells, 1.0);
	p_inv_spacing = 1.0 / (2.2 * p_radius);
	p_num_x = static_cast<int>(std::floor(width * p_inv_spacing)) + 1;
	p_num_y = static_cast<int>(std::floor(height * p_inv_spacing)) + 1;
	int p_num_cells = p_num_x * p_num_y;
	num_cell_particles.resize(p_num_cells, 0);
	first_cell_particle.resize(p_num_cells + 1, 0);
	cell_particle_ids.resize(p_num_cells, 0);
}

double FLIPWaterSimulation::get_cell_size() {
	return cell_size;	
}

double FLIPWaterSimulation::get_particle_radius() {
	return p_radius;	
}

void FLIPWaterSimulation::set_row(const int row, const Array r_cell_types) {
	for (int i = 0; i < f_num_x; i++) {
		double x = i * cell_size;
		double y = row * cell_size;
		int index = i + (f_num_x * row);

		double _s = 1.0;

		int cell_type = r_cell_types[i];

		if (cell_type == SOLID_CELL) {
			cell_type = SOLID_CELL;
			_s = 0.0;
		} else if (cell_type == FLUID_CELL) {
			for (int n = 0; n < 2; n++) {
				double dx = (cell_size / 3) + ((n % 2) * (cell_size / 2));
				double dy = (cell_size / 3) + ((n % 2) * (cell_size / 2));
				p_positions.push_back(Vector2(x + dx, y + dy));
				p_velocities.push_back(Vector2(0.0, 0.0));
			}
		}

		if (row == f_num_y - 1) {
			cell_type = SOLID_CELL;
			_s = 0.0;
		}

		cell_types[index] = cell_type;
		s[index] = _s;
	}
}

void FLIPWaterSimulation::set_cell_types(const Array r_cell_types) {
	for (int i = 0; i < r_cell_types.size(); i++) {
		int cell_type = r_cell_types[i];
		double _s = 1.0;

		if (cell_type == SOLID_CELL) {
			_s = 0.0;
		} else if (cell_type == FLUID_CELL) {
			double x = (i % f_num_x) * cell_size;
			double y = static_cast<int>(std::floor(i / f_num_x)) * cell_size;

			for (int n = 0; n < 2; n++) {
				double dx = (cell_size / 3) + ((n % 2) * (cell_size / 2));
				double dy = (cell_size / 3) + ((n % 2) * (cell_size / 2));
				p_positions.push_back(Vector2(x + dx, y + dy));
				p_velocities.push_back(Vector2(0.0, 0.0));
			}
		}

		cell_types[i] = cell_type;
		s[i] = _s;
	}
}

void FLIPWaterSimulation::add_wall(Vector2 pos) {
	int i = static_cast<int>(std::floor(pos.x / cell_size));
	int j = static_cast<int>(std::floor(pos.y / cell_size));
	int index = i + (f_num_x * j);

	if (index > cell_types.size() - 1) {
		return;	
	}
	
	s[index] = 0.0;
	cell_types[index] = SOLID_CELL;
	p_densities[index] = 0.0;
    pressure[index] = 0.0;
    // u[index] = 0.0;
    // v[index] = 0.0;
	// pressure[index] = 0.0;

	// Reset pressure at the four surrounding corners to remove any residual effects
	int left = (i > 0) ? i - 1 : i;
	int right = (i < f_num_x - 1) ? i + 1 : i;
	int top = (j > 0) ? j - 1 : j;
	int bottom = (j < f_num_y - 1) ? j + 1 : j;

    u[left + (f_num_x * j)] = 0.0;
    u[right + (f_num_x * j)] = 0.0;
    v[i + (f_num_x * top)] = 0.0;
    v[i + (f_num_x * bottom)] = 0.0;

	// Reset pressure at the cell itself and the four adjacent cells
	// pressure[left + (f_num_x * j)] = 0.0;
	// pressure[right + (f_num_x * j)] = 0.0;
	// pressure[i + (f_num_x * top)] = 0.0;
	// pressure[i + (f_num_x * bottom)] = 0.0;

	// Additional step: If you want to ensure neighboring cell influences are minimized,
	// consider adjusting the density field (`p_densities`) in these cells as well:
	// p_densities[index] = 0.0;
	// p_densities[left + (f_num_x * j)] = 0.0;
	//p_densities[right + (f_num_x * j)] = 0.0;
	// p_densities[i + (f_num_x * top)] = 0.0;
	//p_densities[i + (f_num_x * bottom)] = 0.0;

	double min_x = i * cell_size;
	double max_x = (i + 1) * cell_size;
	double min_y = j * cell_size;
	double max_y = (j + 1) * cell_size;

	// Remove or adjust particles within the wall cell bounds
	for (int p = 0; p < p_positions.size(); p++) {
		Vector2 particle_pos = p_positions[p];
		
		if (particle_pos.x >= min_x && particle_pos.x < max_x &&
		    particle_pos.y >= min_y && particle_pos.y < max_y) {

			// Option 1: Remove the particle completely
			// p_positions.erase(p_positions.begin() + p);
			// p_velocities.erase(p_velocities.begin() + p);
			// p--;  // Adjust index after erasing

			// Option 2: Push particle out to the nearest edge of the cell
			// Uncomment the following lines if you prefer repositioning:
			double dx = std::min(particle_pos.x - min_x, max_x - particle_pos.x);
			double dy = std::min(particle_pos.y - min_y, max_y - particle_pos.y);
			if (dx < dy) {
				particle_pos.x += (particle_pos.x - min_x < max_x - particle_pos.x) ? -dx - p_radius : dx + p_radius;
			} else {
				particle_pos.y += (particle_pos.y - min_y < max_y - particle_pos.y) ? -dy - p_radius : dy + p_radius;
			}
			p_positions[p] = particle_pos;
			p_velocities[p] = Vector2(0.0, 0.0);  // Reset velocity
		}
	}
	
	// Calculate the boundaries of the wall cell
    // double min_x = i * cell_size;
    // double max_x = (i + 1) * cell_size;
    // double min_y = j * cell_size;
    // double max_y = (j + 1) * cell_size;

	// std::vector<int> particles_to_remove;

	// // Remove particles in the cell
	// for (int i = 0; i < p_positions.size(); i++) {
    //     if (pos.x >= min_x && pos.x < max_x && pos.y >= min_y && pos.y < max_y) {
	// 		particles_to_remove.insert(particles_to_remove.begin(), i);
	// 	}
	// }

	// if (particles_to_remove.size() > 0) {
	// 	std::sort(particles_to_remove.begin(), particles_to_remove.end(), std::greater<int>());
	// 	for (int i = 0; i < particles_to_remove.size(); i++) {
	// 		int index = particles_to_remove[i];
	// 		p_positions.erase(p_positions.begin() + index);
	// 		p_velocities.erase(p_velocities.begin() + index);
	// 	}
	// }

    // Push particles out of the new wall cell
    // for (auto& pos : p_positions) {
    //     if (pos.x >= min_x && pos.x < max_x && pos.y >= min_y && pos.y < max_y) {
    //         // Find the closest edge of the cell and push the particle to it
    //         double dx = std::min(pos.x - min_x, max_x - pos.x) + 1;
    //         double dy = std::min(pos.y - min_y, max_y - pos.y) + 1;

    //         if (dx < dy) {
    //             pos.x += pos.x - min_x < max_x - pos.x ? -dx : dx;
    //         } else {
    //             pos.y += pos.y - min_y < max_y - pos.y ? -dy : dy;
    //         }

    //         // Optionally, reset the velocity of the moved particle
    //         int index = &pos - &p_positions[0];
    //         p_velocities[index] = Vector2(0.0, 0.0);
    //     }
    // }
}

void FLIPWaterSimulation::remove_wall(Vector2 pos) {
	int i = static_cast<int>(std::floor(pos.x / cell_size));
	int j = static_cast<int>(std::floor(pos.y / cell_size));
	int index = i + (f_num_x * j);
	
	if (i > 0 && i < f_num_x - 1 && index < cell_types.size() - 1) {
		s[index] = 1.0;
		cell_types[index] = AIR_CELL;
	}
}

Array FLIPWaterSimulation::get_particle_data() {
	Array out;
	for (int i = 0; i < p_positions.size(); i++) {
		Dictionary p_data;
		p_data["index"] = i;
		p_data["position"] = p_positions[i];
		p_data["velocity"] = p_velocities[i];
		out.append(p_data);
	}
	return out;
}

Array FLIPWaterSimulation::get_cell_types() {
	Array out;
	for (int i = 0; i < cell_types.size(); i++) {
		out.append(cell_types[i]);
	}
	return out;
}

Array FLIPWaterSimulation::get_pressure() {
	Array out;
	for (int i = 0; i < pressure.size(); i++) {
		out.append(pressure[i]);
	}
	return out;
}

void FLIPWaterSimulation::push_particle(int index, Vector2 vel) {
	p_velocities[index] += vel;
}

void FLIPWaterSimulation::add_particle(Vector2 pos) {
	p_positions.push_back(Vector2(pos.x, pos.y));
	p_velocities.push_back(Vector2(0.0, 0.0));
}

void FLIPWaterSimulation::remove_particle(int index) {
	p_positions.erase(p_positions.begin() + index);
	p_velocities.erase(p_velocities.begin() + index);
}

void FLIPWaterSimulation::add_jet(Vector2 pos, Vector2 dir) {
	jet_positions.push_back(Vector2(pos.x, pos.y));
	jet_directions.push_back(Vector2(dir.x, dir.y));
}

void FLIPWaterSimulation::shift_down() {
	for (int j = f_num_y - 1; j > 0; j--) {
        for (int i = 0; i < f_num_x; i++) {
            int index = i + (f_num_x * j);
            int indexAbove = i + (f_num_x * (j - 1));
			cell_types[index] = cell_types[indexAbove];
			s[index] = s[indexAbove];
			u[index] = u[indexAbove];
			v[index] = v[indexAbove];
			du[index] = du[indexAbove];
			dv[index] = dv[indexAbove];
			prev_u[index] = prev_u[indexAbove];
			prev_v[index] = prev_v[indexAbove];
			pressure[index] = pressure[indexAbove];
			p_densities[index] = p_densities[indexAbove];
        }
    }

    // Handle the new top row (e.g., set to AIR_CELL)
    // for (int i = 0; i < f_num_x; i++) {
    //     int index = i;
    //     cell_types[index] = AIR_CELL;
    //     s[index] = 1.0;
    // }

	for (int i = 0; i < p_positions.size(); i++) {
		p_positions[i].y += cell_size;	
	}

	for (int i = 0; i < jet_positions.size(); i++) {
		jet_positions[i].y += cell_size;	
	}
}

void FLIPWaterSimulation::shift_up() {
	for (int j = 0; j < f_num_y - 1; j++) {
        for (int i = 0; i < f_num_x; i++) {
            int index = i + (f_num_x * j);
            int indexBelow = i + (f_num_x * (j + 1));
			cell_types[index] = cell_types[indexBelow];
			s[index] = s[indexBelow];
			u[index] = u[indexBelow];
			v[index] = v[indexBelow];
			du[index] = du[indexBelow];
			dv[index] = dv[indexBelow];
			prev_u[index] = prev_u[indexBelow];
			prev_v[index] = prev_v[indexBelow];
			pressure[index] = pressure[indexBelow];
			p_densities[index] = p_densities[indexBelow];
        }
    }

    // Handle the new bottom row (e.g., set to AIR_CELL)
    // for (int i = 0; i < f_num_x; i++) {
    //     int index = i + f_num_x * (f_num_y - 1);
    //     cell_types[index] = AIR_CELL;
    //     s[index] = 1.0;
    // }

	for (int i = 0; i < p_positions.size(); i++) {
		p_positions[i].y -= cell_size;	
	}

	for (int i = 0; i < jet_positions.size(); i++) {
		jet_positions[i].y -= cell_size;	
	}
}

int FLIPWaterSimulation::_clamp_int(int v, int min, int max) {
	return std::min(max, std::max(min, v));
}

double FLIPWaterSimulation::_clamp_double(double v, double min, double max) {
	return std::min(max, std::max(min, v));
}

void FLIPWaterSimulation::update(double delta) {
	_integrate_particles(delta);
	_push_particles_apart();
	_handle_collisions();
	_transfer_velocities(true);
	_update_particle_density();
	_solve_incompressibility(delta, 30);
	_transfer_velocities(false);
}

void FLIPWaterSimulation::_integrate_particles(double delta) {
	double max_y = ((f_num_y - 1) * cell_size) - p_radius;
	// std::vector<size_t> indices(p_positions.size());
	// std::iota(indices.begin(), indices.end(), 0);

	// std::transform(std::execution::par, indices.begin(), indices.end(), indices.begin(), [&](size_t i) -> size_t {
    for (int i = 0; i < p_positions.size(); i++) {
		if (p_positions[i].y > max_y) {
			// return i;
            continue;
		}
		p_velocities[i] += Vector2(0.0, 9.81) * delta;

		for (int j = 0; j < jet_positions.size(); j++) {
			Vector2 vector_to_particle = p_positions[i] - jet_positions[j];
            double distance_to_jet = vector_to_particle.length();

            if (distance_to_jet < 1280.0) {
                Vector2 jet_dir = jet_directions[j].normalized();
                Vector2 normalized_vector_to_particle = vector_to_particle.normalized();
                
                double dot_product = jet_dir.dot(normalized_vector_to_particle);
                double angle_radians = std::acos(dot_product);
                
                if (angle_radians < 0.3) { // Assuming the angle threshold is also in radians
                    double force_magnitude = 10.0 * (1.0 - (distance_to_jet / 1280.0));
                    Vector2 force = normalized_vector_to_particle * force_magnitude;
                    p_velocities[i] += force * delta; // Apply force to velocity
                }
            }
		}

		p_velocities[i].x = _clamp_double(p_velocities[i].x, -40.0, 40.0);
		p_velocities[i].y = _clamp_double(p_velocities[i].y, -40.0, 40.0);
		p_positions[i] += p_velocities[i];
		// return i;
    }
	// });
}

void FLIPWaterSimulation::_handle_collisions() {
	double min_x = cell_size + p_radius;
	double max_x = ((f_num_x - 1) * cell_size) - p_radius;
	double min_y = cell_size + p_radius;
	double max_y = ((f_num_y - 1) * cell_size) - p_radius;
	
	std::vector<int> particles_to_remove;
	
	for (int i = 0; i < p_positions.size(); i++) {
		double x = p_positions[i].x;
		double y = p_positions[i].y;
		
		if (x < min_x) {
			// x = min_x;
			// p_velocities[i].x = 0.0;
			particles_to_remove.insert(particles_to_remove.begin(), i);
			continue;
		} else if (x > max_x) {
			// x = max_x;
			// p_velocities[i].x = 0.0;
			particles_to_remove.insert(particles_to_remove.begin(), i);
			continue;
		}
		
		if (y < min_y) {
			particles_to_remove.insert(particles_to_remove.begin(), i);
			continue;
			// p_velocities[i].x = 0.0;
			// p_velocities[i].y = 0.0;
		} else if (y > max_y) {
			// y = min_y;
			// particles_to_remove.insert(particles_to_remove.begin(), i);
			// continue;
			p_velocities[i].x = 0.0;
			p_velocities[i].y = 0.0;
		}
		
		int pi = static_cast<int>(std::floor(x / cell_size));
		int pj = static_cast<int>(std::floor(y / cell_size));
		int index = pi + (f_num_x * pj);
		
		// bool isEnclosed = false;
		// // if (s[index - 1] != 0.0 || s[index + 1] != 0.0 || s[index - f_num_x] != 0.0 || s[index + f_num_x] != 0.0) {
		// if (s[index] == 0.0) {
		// 	// isEnclosed = false;
		// 	isEnclosed = true;
		// }

		// if (isEnclosed) {
		// 	particles_to_remove.insert(particles_to_remove.begin(), i);
		// 	continue; // No further checks needed if the particle is going to be removed
		// }
		
		if (s[index - 1] == 0.0 && x - p_radius < pi * cell_size) {
			x = (pi * cell_size) + p_radius;
			p_velocities[i].x = 0.0;
		} else if (s[index + 1] == 0.0 && x + p_radius > (pi + 1) * cell_size) {
			x = ((pi + 1) * cell_size) - p_radius;
			p_velocities[i].x = 0.0;
		}

		if (s[index - f_num_x] == 0.0 && y - p_radius < pj * cell_size) {
			y = (pj * cell_size) + p_radius;
			p_velocities[i].y = 0.0;
		} else if (s[index + f_num_x] == 0.0 && y + p_radius > (pj + 1) * cell_size) {
			y = ((pj + 1) * cell_size) - p_radius;
			p_velocities[i].y = 0.0;
		}
		
		p_positions[i].x = x;
		p_positions[i].y = y;
	}
	
	if (particles_to_remove.size() > 0) {
		std::sort(particles_to_remove.begin(), particles_to_remove.end(), std::greater<int>());
		for (int i = 0; i < particles_to_remove.size(); i++) {
			int index = particles_to_remove[i];
			p_positions.erase(p_positions.begin() + index);
			p_velocities.erase(p_velocities.begin() + index);
		}
	}
}

void FLIPWaterSimulation::_push_particles_apart() {
	if (p_positions.size() == 0) {
		return;
	}

	std::fill(num_cell_particles.begin(), num_cell_particles.end(), 0);
	int p_num_cells = p_num_x * p_num_y;
	
	for (int i = 0; i < p_positions.size() - 1; i++) {
		Vector2 pos = p_positions[i];
		double x = pos.x;
		double y = pos.y;
		int xi = _clamp_int(static_cast<int>(std::floor(x * p_inv_spacing)), 0, p_num_x - 1);
		int yi = _clamp_int(static_cast<int>(std::floor(y * p_inv_spacing)), 0, p_num_y - 1);
		int cell_nr = xi + (p_num_x * yi);
		num_cell_particles[cell_nr]++;
	}

	int _first = 0;
	for (int i = 0; i < p_num_cells - 1; i++) {
		_first += num_cell_particles[i];
		first_cell_particle[i] = _first;	
	}
	first_cell_particle[p_num_cells] = _first;
	
	for (int i = 0; i < p_positions.size() - 1; i++) {
		Vector2 pos = p_positions[i];
		double x = pos.x;
		double y = pos.y;
		int xi = _clamp_int(static_cast<int>(std::floor(x * p_inv_spacing)), 0, p_num_x - 1);
		int yi = _clamp_int(static_cast<int>(std::floor(y * p_inv_spacing)), 0, p_num_y - 1);
		int cell_nr = xi + (p_num_x * yi);
		first_cell_particle[cell_nr]--;
		cell_particle_ids[first_cell_particle[cell_nr]] = i;
	}
	
	double min_dist = 2.0 * p_radius;
	double min_dist2 = min_dist * min_dist;
	for (int iter = 0; iter < 2; iter++) {
		for (int i = 0; i < p_positions.size() - 1; i++) {
			Vector2 pos = p_positions[i];
			double x = pos.x;
			double y = pos.y;
			
			int pxi = static_cast<int>(std::floor(x * p_inv_spacing));
			int pyi = static_cast<int>(std::floor(y * p_inv_spacing));
			int x0 = std::max(pxi - 1, 0);
			int y0 = std::max(pyi - 1, 0);
			int x1 = std::min(pxi + 1, p_num_x - 1);
			int y1 = std::min(pyi + 1, p_num_y - 1);
			
			for (int xi = x0; xi < x1 + 1; xi++) {
				for (int yi = y0; yi < y1 + 1; yi++) {
					int cell_nr = xi + (p_num_x * yi);
					int first = first_cell_particle[cell_nr];
					int last = first_cell_particle[cell_nr + 1];

					for (int j = first; j < last; j++) {
						int id = cell_particle_ids[j];
						if (id == i) {
							continue;
						}
						double qx = p_positions[id].x;
						double qy = p_positions[id].y;
						
						double dx = qx - x;
						double dy = qy - y;
						double d2 = dx * dx + dy * dy;
						if (d2 > min_dist2 || d2 == 0.0) {
							continue;
						}
						
						double d = std::sqrt(d2);
						double _s = 0.5 * (min_dist - d) / d;
						dx *= _s;
						dy *= _s;
						p_positions[i].x -= dx;
						p_positions[i].y -= dy;
						p_positions[id].x += dx;
						p_positions[id].y += dy;
					}
				}
			}
		}
	}
}

void FLIPWaterSimulation::_transfer_velocities(bool to_grid = false) {
	double max_y = ((f_num_y - 1) * cell_size) - p_radius;
	double h = cell_size;
	double h1 = f_inv_spacing;
	double h2 = 0.5 * h;
	
	if (to_grid == true) {
		prev_u = u;
		prev_v = v;
		std::fill(du.begin(), du.end(), 0.0);
		std::fill(dv.begin(), dv.end(), 0.0);
		std::fill(u.begin(), u.end(), 0.0);
		std::fill(v.begin(), v.end(), 0.0);

		for (int i = 0; i < f_num_cells; i++) {
			cell_types[i] = s[i] == 0.0 ? SOLID_CELL : AIR_CELL;
		}
		
		for (int i = 0; i < p_positions.size(); i++) {
			Vector2 pos = p_positions[i];
			
			if (pos.y > max_y) {
				continue;
			}

			int xi = _clamp_int(static_cast<int>(std::floor(pos.x * h1)), 0, f_num_x - 1);
			int yi = _clamp_int(static_cast<int>(std::floor(pos.y * h1)), 0, f_num_y - 1);
			int cell_nr = xi + (f_num_x * yi);
			if (cell_types[cell_nr] == AIR_CELL) {
				cell_types[cell_nr] = FLUID_CELL;
			}
		}
	}
	
	for (int component = 0; component < 2; component++) {
		double dx = component == 0 ? 0.0 : h2;
		double dy = component == 0 ? h2 : 0.0;

		std::vector<double>& f = component == 0 ? u : v;
		std::vector<double>& prev_f = component == 0 ? prev_u : prev_v;
		std::vector<double>& d = component == 0 ? du : dv;

		for (int i = 0; i < p_positions.size(); i++) {
			Vector2 pos = p_positions[i];
			
			if (pos.y > max_y) {
				continue;
			}
			
			double x = pos.x;
			double y = pos.y;

			x = _clamp_double(x, h, f_num_x * h);
			y = _clamp_double(y, h, f_num_y * h);
			
			int x0 = std::min(static_cast<int>(std::floor((x - dx) * h1)), f_num_x - 1);
			double tx = ((x - dx) - x0 * h) * h1;
			int x1 = std::min(x0 + 1, f_num_x - 1);

			int y0 = std::min(static_cast<int>(std::floor((y - dy) * h1)), f_num_y - 1);
			double ty = ((y - dy) - y0 * h) * h1;
			int y1 = std::min(y0 + 1, f_num_y - 1);
			
			double sx = 1.0 - tx;
			double sy = 1.0 - ty;
			
			double d0 = sx * sy;
			double d1 = tx * sy;
			double d2 = tx * ty;
			double d3 = sx * ty;

			int nr0 = x0 + (f_num_x * y0);
			int nr1 = x1 + (f_num_x * y0);
			int nr2 = x1 + (f_num_x * y1);
			int nr3 = x0 + (f_num_x * y1);
			
			double v = component == 0 ? p_velocities[i].x : p_velocities[i].y;
			
			if (to_grid == true) {
				if (cell_types[nr0] != SOLID_CELL) {
					f[nr0] += v * d0;
					d[nr0] += d0;
				}
				if (cell_types[nr1] != SOLID_CELL) {
					f[nr1] += v * d1;
					d[nr1] += d1;
				}
				if (cell_types[nr2] != SOLID_CELL) {
					f[nr2] += v * d2;
					d[nr2] += d2;
				}
				if (cell_types[nr3] != SOLID_CELL) {
					f[nr3] += v * d3;
					d[nr3] += d3;
				}
			} else {
				int offset = component == 1 ? f_num_x : 1;
				double valid0 = cell_types[nr0] != SOLID_CELL && cell_types[nr0 - offset] != SOLID_CELL ? 1.0 : 0.0;
				double valid1 = cell_types[nr1] != SOLID_CELL && cell_types[nr1 - offset] != SOLID_CELL ? 1.0 : 0.0;
				double valid2 = cell_types[nr2] != SOLID_CELL && cell_types[nr2 - offset] != SOLID_CELL ? 1.0 : 0.0;
				double valid3 = cell_types[nr3] != SOLID_CELL && cell_types[nr3 - offset] != SOLID_CELL ? 1.0 : 0.0;
				double _d = (valid0 * d0) + (valid1 * d1) + (valid2 * d2) + (valid3 * d3);
				
				double flip_ratio = 0.9;

				if (_d > 0.0) {
					double pic_v = ((valid0 * d0 * f[nr0]) + (valid1 * d1 * f[nr1]) + (valid2 * d2 * f[nr2]) + (valid3 * d3 * f[nr3])) / _d;
					double corr = ((valid0 * d0 * (f[nr0] - prev_f[nr0])) + (valid1 * d1 * (f[nr1] - prev_f[nr1])) + (valid2 * d2 * (f[nr2] - prev_f[nr2])) + (valid3 * d3 * (f[nr3] - prev_f[nr3]))) / _d;
					double flip_v = v + corr;

					if (component == 0) {
						p_velocities[i].x = ((1.0 - flip_ratio) * pic_v) + flip_ratio * flip_v;
					} else {
						p_velocities[i].y = ((1.0 - flip_ratio) * pic_v) + flip_ratio * flip_v;
					}
				}
			}
		}
		
		if (to_grid == true) {
			for (int i = 0; i < f.size(); i++) {
				if (d[i] > 0.0) {
					f[i] /= d[i];
				}
			}
			
			for (int i = 0; i < f_num_x; i++) {
				for (int j = 0; j < f_num_y; j++) {
					bool solid = cell_types[i + (f_num_x * j)] == SOLID_CELL;
					if (solid || (i > 0 && cell_types[i - 1 + (f_num_x * j)] == SOLID_CELL)) {
						u[i + (f_num_x * j)] = prev_u[i + (f_num_x * j)];
					}
					if (solid || (j > 0 && cell_types[i + (f_num_x * (j - 1))] == SOLID_CELL)) {
						v[i + (f_num_x * j)] = prev_v[i + (f_num_x * j)];
					}
				}
			}
		}
	}
}

void FLIPWaterSimulation::_update_particle_density() {
	double max_y = ((f_num_y - 1) * cell_size) - p_radius;
	std::fill(p_densities.begin(), p_densities.end(), 0.0);
	
	for (int i = 0; i < p_positions.size(); i++) {
		Vector2 pos = p_positions[i];
		
		if (pos.y > max_y) {
			continue;
		}

		double x = pos.x;
		double y = pos.y;
		
		x = _clamp_double(x, cell_size, (f_num_x - 1) * cell_size);
		y = _clamp_double(y, cell_size, (f_num_y - 1) * cell_size);
		
		int x0 = static_cast<int>(std::floor((x - (0.5 * cell_size)) * f_inv_spacing));
		double tx = ((x - (0.5 * cell_size)) - x0 * cell_size) * f_inv_spacing;
		int x1 = std::min(x0 + 1, f_num_x - 1);

		int y0 = static_cast<int>(std::floor((y - (0.5 * cell_size)) * f_inv_spacing));
		double ty = ((y - (0.5 * cell_size)) - y0 * cell_size) * f_inv_spacing;
		int y1 = std::min(y0 + 1, f_num_y - 1);
		
		double sx = 1.0 - tx;
		double sy = 1.0 - ty;
		
		if (x0 < f_num_x && y0 < f_num_y) {
			p_densities[x0 + (f_num_x * y0)] += sx * sy;
		}
		if (x1 < f_num_x && y0 < f_num_y) {
			p_densities[x1 + (f_num_x * y0)] += tx * sy;
		}
		if (x1 < f_num_x && y1 < f_num_y) {
			p_densities[x1 + (f_num_x * y1)] += tx * ty;
		}
		if (x0 < f_num_x && y1 < f_num_y) {
			p_densities[x0 + (f_num_x * y1)] += sx * ty;
		}
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

void FLIPWaterSimulation::_solve_incompressibility(double delta, int iters) {
	prev_u = u;
	prev_v = v;

	double cp = density * cell_size / delta;

	for (int iter = 0; iter < iters; iter++) {
		for (int i = 1; i < f_num_x - 1; i++) {
			for (int j = 1; j < f_num_y - 1; j++) {
				if (cell_types[i + (f_num_x * j)] != FLUID_CELL) {
					continue;
				}
				
				int center = i + (f_num_x * j);
				int left = i + (f_num_x * j) - 1;
				int right = i + (f_num_x * j) + 1;
				int top = i + (f_num_x * (j - 1));
				int bottom = i + (f_num_x * (j + 1));

				double s0 = s[center];
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
					double k = 1.0;
					double compression = p_densities[i + (f_num_x * j)] - p_rest_density;
					if (compression > 0.0) {
						div = div - (k * compression);
					}
				}
				
				double _p = -div / _s;
				_p *= 1.9;
				pressure[center] += cp * _p;

				// if (s[left] == 0.0 || s[right] == 0.0 || s[top] == 0.0 || s[bottom] == 0.0) {
				// 	_p *= 0.5;  // Reduce effect if wall is nearby
				// }

				u[center] -= sx0 * _p;
				u[right] += sx1 * _p;
				v[center] -= sy0 * _p;
				v[bottom] += sy1 * _p;
			}
		}
	}	
}