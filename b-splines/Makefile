JULIA="/home/mn/Developer/working-copies/julia/julia"

all: drawings

drawings:
	$(JULIA) --eval "include(\"bspline.jl\"); \
		exercise_zero(); \
		exercise_six();	";
	gnuplot bspline.gnuplot;
