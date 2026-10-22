# FLIP Water Simulation (Godot 4)

A 2D FLIP/PIC fluid simulation for Godot 4.2+, written as a GDExtension in C++. You drop in a `FLIPWaterSimulation` node and it collides with your existing physics bodies and TileMaps, renders itself, and can be driven by emitter, drain, force and sensor nodes.

## Building

```sh
git submodule update --init
scons                      # debug build for your platform
scons target=template_release
```

You need SCons and a C++17 compiler. The library is written to `demo/addons/flip_water_simulation/bin/`.

## Installing

Copy `demo/addons/flip_water_simulation/` (including the built `bin/`) into your project's `addons/` folder. Godot loads GDExtensions automatically, so there's nothing to enable under Project Settings → Plugins.

`demo/` is a runnable example. Left drag moves the paddle, right drag paints walls, middle click spawns water, and Space cycles the render modes.

## Quick start

1. Add a `FLIPWaterSimulation` node. Its domain is `grid_size × cell_size` pixels, starting at the node's origin. The editor shows the domain as an outline.
2. Add water from code:
   ```gdscript
   $Water.fill_cells(Rect2i(2, 20, 16, 15), FLIPWaterSimulation.FLUID_CELL)
   ```
   Alternatively, add a `FLIPEmitter` child.
3. Run. Any `StaticBody2D`, `AnimatableBody2D`, `CharacterBody2D`, `RigidBody2D` or `TileMap` on the `collision_mask` layers acts as a wall. Moving bodies push the water.

## Nodes

### FLIPWaterSimulation

| Group | Properties |
| --- | --- |
| Grid | `grid_size` (cells), `cell_size` (px), `edge_left` / `edge_top` / `edge_right` / `edge_bottom`: `Wall` keeps water in, `Open` lets it leave and deletes it |
| Simulation | `gravity` (px/s², default `(0, 980)`), `flip_ratio` (0 = viscous PIC, 1 = splashy FLIP), `max_speed` (0 = only the one-cell-per-step limit), `solver_iterations`, `substeps`, `max_particles`, `auto_update`, `use_threads` |
| Scrolling | `follow_target`: the domain re-centres on this node in whole-cell steps |
| Collision | `collision_mask`: physics layers that act as walls |
| Rendering | `render_mode` (`None`, `Particles`, `Surface`), `water_color`, `foam_color`, `foam_speed` (speed at which water turns fully foam-coloured), `particle_scale`, `surface_resolution` |
| Debug | `debug_draw` shows solid, moving and fluid cells |

| Method | Description |
| --- | --- |
| `set_cell(cell, type)`, `fill_cells(rect, type)`, `get_cell(cell)` | Paint `SOLID_CELL` / `AIR_CELL`, or spawn water with `FLUID_CELL`. |
| `local_to_cell(pos)`, `cell_to_local(cell)`, `get_domain_rect()` | Coordinate helpers (node-local space). |
| `shift(cells)` | Move the domain by whole cells. The water stays put in world space, and water that ends up outside the domain is deleted. Emits `region_exposed(rect)` for each newly uncovered strip. |
| `add_particle(pos, velocity)`, `apply_impulse(center, radius, velocity)`, `remove_particles_in_circle(center, radius)` | Add, push or remove water. |
| `count_particles_in_rect(rect)`, `count_particles_in_circle(center, radius)` | Gameplay queries. |
| `get_particle_positions()`, `get_particle_velocities()`, `get_cell_types()`, `get_pressure()` | Raw state as packed arrays, e.g. for custom rendering. |
| `update(delta)` | Step manually (set `auto_update = false`). |
| `refresh_colliders()` | Re-read static bodies (see below). |
| `clear()` | Remove all water and painted walls. |

### Effectors

Effectors act on the simulation set in their `simulation` property. If that's empty, they use the nearest `FLIPWaterSimulation` ancestor. Each has an `enabled` toggle and draws a gizmo in the editor.

| Node | Properties | Notes |
| --- | --- | --- |
| `FLIPEmitter` | `rate` (/s), `speed`, `spread`, `width` | Emits along its +X axis. `emit_burst(count)` spawns extra water on the next step. |
| `FLIPDrain` | `radius` | Removes water. Emits `water_drained(count)`. |
| `FLIPForce` | `mode` (`MODE_RADIAL` / `MODE_DIRECTIONAL`), `strength` (px/s², negative attracts), `radius`, `spread`, `falloff` | Directional mode pushes along +X, and `spread` limits the effect to a cone. |
| `FLIPSensor` | `size`, `threshold` | `get_particle_count()` and `is_wet()`. Emits `water_entered` and `water_exited`. |

## Notes

- **Colliders** are rasterised onto the grid every frame, so precision is one cell, and shapes thinner than about half a cell can leak. Static bodies are cached for speed. New ones are picked up automatically, but call `refresh_colliders()` after moving or removing a `StaticBody2D`.
- **Scrolling:** with `follow_target`, use `Open` edges and place effectors outside the simulation node, otherwise they scroll along with it.
- **Performance:** the solver uses Godot's `WorkerThreadPool` for large grids and particle counts, and gives identical results with threading on or off. A 240×140 grid with 26k particles takes about 4 ms per frame on a 16-thread CPU, versus 7 ms on one thread.
- **Surface mode** draws particles into an internal `SubViewport` and thresholds the result. Its cost scales with `surface_resolution`.
