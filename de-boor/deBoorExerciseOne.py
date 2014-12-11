

from deBoorDrawing import * 

def exercise_one():
    """
    Simple exercise to plot an open mushroom.
    """

    interval = (0,1)

    control_net = np.matrix([
                             [-0.2, 2],
                             [-0.3, 6.2],
                             [-1.2, 4.8],
                             [-2.8, 8.8],
                             [-0.7, 14],
                             [1.4, 14.7],
                             [3.6, 10.2],
                             [3.2, 5.1],
                             [1.5, 6.2],
                             [1.4, 2],
                             ])

#   First we plot a curve where internal knots have maximum multiplicities
    arguments = {   
        'order':4, 
        'interval':interval,
        'internal_knots':sample_internal_knots_uniformly_in(interval, 3),
        'control_net':control_net,
        'multiplicities':[2,2,2]
    }

    curve = draw(**arguments)


#   After we plot a curve where each internal knot have multiplicity 1.
    plot_curve(curve, control_net, axis=[-4, 4, 0, 16])

    arguments = {   
        'order':4, 
        'interval':interval,
        'internal_knots':sample_internal_knots_uniformly_in(interval, 6),
        'control_net':control_net,
    }

    curve = draw(**arguments)

    plot_curve(curve, control_net, axis=[-4, 4, 0, 16])


#________________________________________________________________________
exercise_one() # run the exercise


