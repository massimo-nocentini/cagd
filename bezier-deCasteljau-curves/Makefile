JULIA="/home/mn/Developer/working-copies/julia/julia"

all: drawings

drawings:
	$(JULIA) --eval "include(\"deCasteljau.jl\"); \
		include(\"exercises.jl\"); \
		exercise_one(); \
		exercise_two(); \
		exercise_four(); \
		exercise_five(); \
		exercise_six(); \
		exercise_seven(); \
		exercise_seven_left(); ";
	gnuplot bezier.gnuplot;
