# JULIA="/home/mn/Developer/working-copies/julia/julia"
JULIA=julia

all: drawings

drawings:
	$(JULIA) --eval "include(\"bspline.jl\"); \
		include(\"exercises.jl\"); \
		exercise_zero(); \
		exercise_six();	\
		exercise_one(); \
		exercise_three(); \
		exercise_four(); \
		exercise_five(); \
		exercise_two(); \
		exercise_knotsInsertion(); ";
	gnuplot bspline.gnuplot;
