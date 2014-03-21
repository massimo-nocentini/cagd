set terminal postscript eps enhanced color font 'Helvetica,10'

set output 'exercise-one.eps'
plot [0:12] [0:10] "bezier-exercise-one.coordinates" with lines, "simple-control-polygon.dat" with linespoints pt 6

set output 'exercise-two.eps'
plot [0.5:3.5] [-0.3:1.3] "bezier-exercise-two.coordinates" with lines, "control-poly-exercise-two.coordinates" with linespoints pt 6

set output 'exercise-four.eps'
plot [-0.5:6.5] [-0.3:5] "control-poly-exercise-four.coordinates" with linespoints pt 6, \
	"cyan-control-poly-exercise-four.coordinates" with linespoints pt 6, \
	"magenta-control-poly-exercise-four.coordinates" with linespoints pt 6, \
	"bezier-diagonal-exercise-four.coordinates" with lines, \
	"bezier-last-exercise-four.coordinates" with lines


set output 'exercise-six-original.eps'
plot "bezier-original-exercise-six.coordinates" with lines, "control-poly-exercise-six.coordinates" with linespoints pt 6

set output 'exercise-six-one-more-degree.eps'
plot "bezier-one-more-degree-exercise-six.coordinates" with lines, "control-poly-exercise-six-one-more-degree.coordinates" with linespoints pt 6

set output 'exercise-six-one-more-degree-comparison.eps'
plot [0:11] [0:7] "bezier-one-more-degree-exercise-six.coordinates" with lines, \
	"bezier-two-more-degree-exercise-six.coordinates" with lines, \
	"bezier-three-more-degree-exercise-six.coordinates" with lines, \
	"bezier-original-exercise-six.coordinates" with lines

set output 'exercise-six-higher-degree-control-poly.eps'
plot [0:11] [0:7] "control-poly-exercise-six.coordinates" with lines, \
	"control-poly-exercise-six-one-more-degree.coordinates" with lines, \
	"control-poly-exercise-six-two-more-degree.coordinates" with lines, \
	"control-poly-exercise-six-three-more-degree.coordinates" with lines, \
	"bezier-original-exercise-six.coordinates" with lines

set output 'exercise-five.eps'
plot "bezier-original-exercise-five.coordinates" with lines, \
	"bezier-repeated-exercise-five.coordinates" with lines

