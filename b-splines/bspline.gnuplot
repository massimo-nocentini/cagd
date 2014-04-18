set terminal postscript eps enhanced color font 'Helvetica,10'

set output 'exercise-zero-clumped.eps'
plot [-3:4]  "exercise-zero-clumped-control-poly.coordinates" with linespoints pt 6, \
     "exercise-zero-clumped-bspline.coordinates" with lines

set output 'exercise-zero-uniformed.eps'
plot [-3:4] "exercise-zero-uniformed-control-poly.coordinates" with linespoints pt 6, \
     "exercise-zero-uniformed-bspline.coordinates" with lines

set output 'exercise-zero-closed.eps'
plot [-3:4] [1:16]"exercise-zero-closed-control-poly.coordinates" with linespoints pt 6, \
     "exercise-zero-closed-bspline.coordinates" with lines

set output 'exercise-six-first-closed.eps'
plot [-3:4] [2:16]"exercise-six-first-closed-control-poly.coordinates" with linespoints pt 6, \
     "exercise-six-first-closed-bspline.coordinates" with lines

set output 'exercise-six-second-closed.eps'
plot [-3:4] [2:16]"exercise-six-second-closed-control-poly.coordinates" with linespoints pt 6, \
     "exercise-six-second-closed-bspline.coordinates" with lines

set output 'exercise-one-clumped.eps'
plot [-4:5] [-3:3] "exercise-one-control-poly.coordinates" with linespoints pt 6, \
     "exercise-one-clumped-bspline-k-two.coordinates" with lines, \
     "exercise-one-clumped-bspline-k-three.coordinates" with lines, \
     "exercise-one-clumped-bspline-k-four.coordinates" with lines, \
     "exercise-one-clumped-bspline-k-five.coordinates" with lines, \
     "exercise-one-clumped-bspline-k-six.coordinates" with lines               

set output 'exercise-three.eps'
plot [-4:6] [-3:4]"exercise-three-control-poly.coordinates" with linespoints pt 6, \
     "exercise-three-bspline.coordinates" with lines

set output 'exercise-four.eps'
plot [-4:4] [1:7]"exercise-four-control-poly.coordinates" with linespoints pt 6, \
     "exercise-four-bspline-first-partition.coordinates" with lines, \
     "exercise-four-bspline-second-partition.coordinates" with lines     

set output 'exercise-five.eps'
plot [-3:3] [-2:2]"exercise-five-control-poly.coordinates" with points pt 6, \
     "exercise-five-bspline-single.coordinates" with lines, \
     "exercise-five-bspline-doubled.coordinates" with lines, \
     "exercise-five-bspline-tripled.coordinates" with lines          

set output 'exercise-two.eps'
set multiplot layout 5,1
plot [-2:3] [-1.5:2]"exercise-two-control-poly.coordinates" with linespoints pt 6, \
     "exercise-two-bspline-first.coordinates" with lines
plot [-2:3] [-1.5:2]"exercise-two-control-poly.coordinates" with linespoints pt 6, \
     "exercise-two-bspline-second.coordinates" with lines
plot [-2:3] [-1.5:2]"exercise-two-control-poly.coordinates" with linespoints pt 6, \
     "exercise-two-bspline-three.coordinates" with lines
plot [-2:3] [-1.5:2]"exercise-two-control-poly.coordinates" with linespoints pt 6, \
     "exercise-two-bspline-four.coordinates" with lines
plot [-2:3] [-1.5:2]"exercise-two-control-poly.coordinates" with linespoints pt 6, \
     "exercise-two-bspline-five.coordinates" with lines
unset multiplot

set output 'exercise-knots-insertion.eps'
plot [-1:5] [-1:3] "exercise-knots-insertion.coordinates" with lines, \
     "exercise-knots-insertion-original-control-points.coordinates" with linespoints pt 6
