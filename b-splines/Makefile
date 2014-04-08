JULIA="/home/mn/Developer/working-copies/julia/julia"

all: drawings

drawings:
	$(JULIA) --eval "include(\"bspline.jl\"); \
		exercise_zero(); \
		exercise_six();	\
		exercise_one(); \
		exercise_three(); \
		exercise_four(); \
		exercise_five(); \
		exercise_two(); ";
	gnuplot bspline.gnuplot;
