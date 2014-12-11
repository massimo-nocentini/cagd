
from deBoorDrawing import * 

def exercise_two():
    """
    This is a simple exercise to plot a closed mushroom.
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
                             [1.4, 2]
                             ])

    arguments = {   
        'order':4, 
        'interval':interval,
        'internal_knots':sample_internal_knots_uniformly_in(interval, 9),
        'control_net':control_net,
        'closed':True
    }

    curve = draw(**arguments)

    plot_curve(curve, close_control_net(control_net), axis=[-4, 4, 0, 16])

#________________________________________________________________________
exercise_two()
