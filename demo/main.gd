extends Node2D

const PADDLE_SPEED := 1200.0

@onready var water: FLIPWaterSimulation = $Water
@onready var paddle: AnimatableBody2D = $Paddle
@onready var sensor: FLIPSensor = $Water/Sensor
@onready var basin_label: Label = $BasinLabel

var walls := {}


func _ready() -> void:
	water.fill_cells(Rect2i(2, 20, 16, 15), FLIPWaterSimulation.FLUID_CELL)
	sensor.water_entered.connect(func(): basin_label.modulate = Color.AQUAMARINE)
	sensor.water_exited.connect(func(): basin_label.modulate = Color.WHITE)


func _physics_process(delta: float) -> void:
	var mouse := get_global_mouse_position()
	if Input.is_mouse_button_pressed(MOUSE_BUTTON_LEFT):
		paddle.global_position = paddle.global_position.move_toward(mouse, PADDLE_SPEED * delta)
	if Input.is_mouse_button_pressed(MOUSE_BUTTON_RIGHT):
		var cell := water.local_to_cell(water.to_local(mouse))
		if not walls.has(cell):
			walls[cell] = true
			water.set_cell(cell, FLIPWaterSimulation.SOLID_CELL)
			queue_redraw()
	if Input.is_mouse_button_pressed(MOUSE_BUTTON_MIDDLE):
		for i in 4:
			water.add_particle(water.to_local(mouse) + Vector2(randf_range(-8, 8), randf_range(-8, 8)))
	basin_label.text = "Basin: %d particles" % sensor.get_particle_count()


func _unhandled_input(event: InputEvent) -> void:
	if event is InputEventKey and event.pressed and event.keycode == KEY_SPACE:
		water.render_mode = (water.render_mode + 1) % 3


func _draw() -> void:
	var size := water.cell_size
	for cell in walls:
		draw_rect(Rect2(Vector2(cell) * size, Vector2(size, size)), Color(0.45, 0.45, 0.5))
