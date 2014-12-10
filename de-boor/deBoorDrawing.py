
from deBoorCore import * 

def draw(order, interval, internal_knots, control_net, 
            tabs=None, points=1000, closed=False, 
            multiplicities=None, extended_vector=None):
    """
    Produces a set of point representing the BSpline interpolation of the given control net.

    In particular, this is a decorator of `deBoor` function since it carries over
    some setup stuff whether the desired curve is closed or not, handling the
    addition of necessary nodes in a smart way. Moreover it checks if the
    given control net, order and multiplicities arguments are coherent
    in order to build the requested curve.
    
    If argument `multiplicities` is None, then assign to each knot belonging
    to internal partition a multiplicity of 1.

    If argument `tabs` is None, then take a uniformly spaced tabulation array.

    If argument `extended_vector` is None, just extend the given internal
    knots partition using `extend_knots_vector` function.
    """

    if multiplicities is None:
        multiplicities = np.repeat(1, len(internal_knots))

    n, d = np.shape(control_net)
    if closed:
        assert n == order + sum(multiplicities) - order + 1
        control_net = np.concatenate((control_net, control_net[:order-1, :]), axis=0)
    else:
        assert n == order + sum(multiplicities)

    if tabs is None:
        a, b = interval
        tabs = np.linspace(start=a, stop=b, num=points*(b-a))

    if extended_vector is None:
        extended_vector = extend_knots_vector(
            order, interval, internal_knots, closed, multiplicities)

    return de_Boor(extended_vector, control_net, tabs)


def plot_curve(curve, control_net=None, axis="image"):
    """
    Plot the given curve, against its control net and saving it if desired.
    """

    import matplotlib.pyplot as plt
    
    if control_net is not None:
        plt.plot(control_net[:,0], control_net[:,1], "o--")

    X, Y = curve[:,0], curve[:,1]
    plt.plot(X, Y)

    plt.axis(axis)
    plt.show()
    # saving the image for now will wait


def close_control_net(control_net, axis=0):
    return np.concatenate((control_net, control_net[0, :]), axis=axis)


def sample_internal_knots_uniformly_in(interval, number_of_knots):
    """
    Returns an numpy array containing internal knots uniformly placed.

    This function is useful to generate a partion of `internal` knots
    to be paired with desired multiplicities. It is implemented to be
    used by `client` functions written by users.

    Simple example sampling from (0,1) segment:
    >>> sample_internal_knots_uniformly_in(interval=(0,1), number_of_knots=3)
    array([ 0.25,  0.5 ,  0.75])

    Another simple test:
    >>> sample_internal_knots_uniformly_in(interval=(0,11), number_of_knots=10)
    array([  1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.])
    """

    a, b = interval
    interval_with_extrema = np.linspace(
        start=a, stop=b, num=number_of_knots+2, endpoint=True)
    return interval_with_extrema[1:-1] # discard the first and the last



