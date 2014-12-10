
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

    arguments = {   
        'order':4, 
        'interval':interval,
        'internal_knots':sample_internal_knots_uniformly_in(interval, 3),
        'control_net':control_net,
        'multiplicities':[2,2,2]
    }

    curve = draw(**arguments)

    plot_curve(curve, control_net, axis=[-4, 4, 0, 16])



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


def exercise_three():
    """
    Insert knot till max smoothness for an open mushroom.
    """
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

    interval = (0,1)
    internal_knots = sample_internal_knots_uniformly_in(interval, 6)
    multiplicities = np.ones(6)
    order, closed, axis = 4, False, [-4, 4, 0, 16]
    
    extended_vector = extend_knots_vector(
        order, interval, internal_knots, closed, multiplicities)

    arguments = {
        "order":order,
        "interval":interval,
        "internal_knots":internal_knots,
        "closed":closed,
        "control_net":control_net,
        "multiplicities":multiplicities,
        "extended_vector":extended_vector
    }

#   First build the raw curve where each knots has multiplicity 1
    curve = draw(**arguments)
    plot_curve(curve, control_net, axis=axis)

    def plot_raising_step(multiplicities, extended_vector, control_net):
        arguments = {
            "order":order,
            "interval":interval,
            "internal_knots":internal_knots,
            "closed":closed,
            "control_net":control_net,
            "multiplicities":multiplicities,
            "extended_vector":extended_vector
        }
        curve = draw(**arguments)
        plot_curve(curve, control_net, axis=axis)
        
    raise_internal_knots_to_max_smooth(
        order, internal_knots, multiplicities, extended_vector, 
        control_net, on_each_rising=plot_raising_step)


def exercise_four():
    """
    Insert knot till max smoothness for an closed mushroom.
    """
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

    interval = (0,1)
    internal_knots = sample_internal_knots_uniformly_in(interval, 9)
    multiplicities = np.ones(9)
    order, closed, axis = 4, True, [-4, 4, 0, 16]
    
    extended_vector = extend_knots_vector(
        order, interval, internal_knots, closed, multiplicities)

    arguments = {
        "order":order,
        "interval":interval,
        "internal_knots":internal_knots,
        "closed":closed,
        "control_net":control_net,
        "multiplicities":multiplicities,
        "extended_vector":extended_vector
    }

#   First build the raw curve where each knots has multiplicity 1
    curve = draw(**arguments)
    plot_curve(curve, control_net, axis=axis)

    def plot_raising_step(multiplicities, extended_vector, control_net):
        arguments = {
            "order":order,
            "interval":interval,
            "internal_knots":internal_knots,
            "closed":closed,
            "control_net":control_net,
            "multiplicities":multiplicities,
            "extended_vector":extended_vector
        }
        curve = draw(**arguments)
        plot_curve(curve, control_net, axis=axis)
        
    raise_internal_knots_to_max_smooth(
        order, internal_knots, multiplicities, extended_vector, 
        control_net, on_each_rising=plot_raising_step)
