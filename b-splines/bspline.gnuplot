set terminal postscript eps enhanced color font 'Helvetica,10'

set output 'exercise-zero-clamped.eps'
plot [-3:4]  "exercise-zero-clamped-control-poly.coordinates" with linespoints pt 6, \
     "exercise-zero-clamped-bspline.coordinates" with lines

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



# set output 'exercise-two.eps'
# plot [0.5:3.5] [-0.3:1.3] "bezier-exercise-two.coordinates" with lines, "control-poly-exercise-two.coordinates" with linespoints pt 6

# set output 'exercise-two-three-axes.eps'
# splot  "bezier-three-axes-exercise-two.coordinates" with lines, "control-three-axes-poly-exercise-two.coordinates" with linespoints pt 6

# set output 'exercise-four.eps'
# plot [-0.5:6.5] [-0.3:5] "control-poly-exercise-four.coordinates" with linespoints pt 6, \
# 	"cyan-control-poly-exercise-four.coordinates" with linespoints pt 6, \
# 	"magenta-control-poly-exercise-four.coordinates" with linespoints pt 6, \
# 	"bezier-diagonal-exercise-four.coordinates" with lines, \
# 	"bezier-last-exercise-four.coordinates" with lines


# set output 'exercise-six-original.eps'
# plot "bezier-original-exercise-six.coordinates" with lines, "control-poly-exercise-six.coordinates" with linespoints pt 6

# set output 'exercise-six-one-more-degree.eps'
# plot "bezier-one-more-degree-exercise-six.coordinates" with lines, "control-poly-exercise-six-one-more-degree.coordinates" with linespoints pt 6

# set output 'exercise-six-one-more-degree-comparison.eps'
# plot [0:11] [0:7] "bezier-one-more-degree-exercise-six.coordinates" with lines, \
# 	"bezier-two-more-degree-exercise-six.coordinates" with lines, \
# 	"bezier-three-more-degree-exercise-six.coordinates" with lines, \
# 	"bezier-original-exercise-six.coordinates" with lines

# set output 'exercise-six-higher-degree-control-poly.eps'
# plot [0:11] [0:7] "control-poly-exercise-six.coordinates" with lines, \
# 	"control-poly-exercise-six-one-more-degree.coordinates" with lines, \
# 	"control-poly-exercise-six-two-more-degree.coordinates" with lines, \
# 	"control-poly-exercise-six-three-more-degree.coordinates" with lines, \
# 	"bezier-original-exercise-six.coordinates" with lines

# set output 'exercise-five.eps'
# plot "bezier-original-exercise-five.coordinates" with lines, \
# 	"bezier-repeated-exercise-five.coordinates" with lines

# set output 'exercise-seven-continuity.eps'
# plot [-0.5:25] [-10.5:4] "control-poly-base-exercise-seven-exercise-seven.coordinates" with linespoints pt 6, \
# 	"control-poly-continuity-exercise-seven.coordinates" with linespoints pt 6, \
# 	"bezier-base-exercise-seven.coordinates" with lines, \
# 	"bezier-continuity-exercise-seven.coordinates" with lines

# set output 'exercise-seven-tangent.eps'
# plot [-0.5:25] [-10.5:4] "control-poly-base-exercise-seven-exercise-seven.coordinates" with linespoints pt 6, \
# 	"control-poly-tangent-exercise-seven.coordinates" with linespoints pt 6, \
# 	"bezier-base-exercise-seven.coordinates" with lines, \
# 	"bezier-tangent-exercise-seven.coordinates" with lines

# set output 'exercise-seven-obsculating.eps'
# plot [-0.5:25] [-10.5:4] "control-poly-base-exercise-seven-exercise-seven.coordinates" with linespoints pt 6, \
# 	"control-poly-obsculating-exercise-seven.coordinates" with linespoints pt 6, \
# 	"bezier-base-exercise-seven.coordinates" with lines, \
# 	"bezier-obsculating-exercise-seven.coordinates" with lines

# set output 'exercise-seven-a_succ_i.eps'
# plot "control-poly-base-obsculating-exercise-seven.coordinates" with linespoints pt 6, \
# 	"control-poly-base-obsculating-with-a-succ-i-exercise-seven.coordinates" with linespoints pt 6

# set output 'exercise-seven-continuity-left.eps'
# plot [-0.5:25] [-5:20] "control-poly-base-left-exercise-seven-exercise-seven.coordinates" with linespoints pt 6, \
# 	"control-poly-continuity-left-exercise-seven.coordinates" with linespoints pt 6, \
# 	"bezier-base-left-exercise-seven.coordinates" with lines, \
# 	"bezier-continuity-left-exercise-seven.coordinates" with lines

# set output 'exercise-seven-tangent-left.eps'
# plot [-0.5:25] [-5:20] "control-poly-base-left-exercise-seven-exercise-seven.coordinates" with linespoints pt 6, \
# 	"control-poly-tangent-left-exercise-seven.coordinates" with linespoints pt 6, \
# 	"bezier-base-left-exercise-seven.coordinates" with lines, \
# 	"bezier-tangent-left-exercise-seven.coordinates" with lines

# set output 'exercise-seven-obsculating-left.eps'
# plot [-0.5:25] [-5:20] "control-poly-base-left-exercise-seven-exercise-seven.coordinates" with linespoints pt 6, \
# 	"control-poly-obsculating-left-exercise-seven.coordinates" with linespoints pt 6, \
# 	"bezier-base-left-exercise-seven.coordinates" with lines, \
# 	"bezier-obsculating-left-exercise-seven.coordinates" with lines

# set output 'exercise-seven-a_i-left.eps'
# plot "control-poly-base-obsculating-left-exercise-seven.coordinates" with linespoints pt 6, \
# 	"control-poly-base-obsculating-with-a-i-left-exercise-seven.coordinates" with linespoints pt 6
