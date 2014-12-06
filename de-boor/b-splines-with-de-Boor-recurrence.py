

import numpy as np

def extend_knots_vector(order, interval, internal_knots, closed=False, multiplicities=None):
    """Produces an extended knots vector for BSpline functions given a simple one.

    This function extend the given knots partition adding auxiliary knots 
    on the left and on the right of the given ones in order to be the support
    for drawing curves, both with `open` control nets both with `closed` ones.
    Fails if the given knots partition doesn't match their multiplicities 
    respect length of their lists.

    Parameters
    ==========
    order -- the order of the curve which will be drawn on this extended knots partition 
    interval -- interval containing the given knot partition. Pay attention that
                interval extrema shouldn't appear in the knots partition. It should
                be given in binary tuple form, ie. (a, b) for some a and b.
    internal_knots -- knots partition: this should be thought as a set, therefore
                        knot multiplicity is specified in an auxiliary vector
    closed --   boolean value that indicate if this knots partition will be used
                to draw closed or open curves. Defaults to False.
    multiplicities --   list of multiplicities, one for each knot specified in 
                        `internal_knots` argument. Defaults to a list of 1s.

    Examples
    ========

    Simple extension for an open curve:
    >>> extend_knots_vector(4, (3, 6), [4, 5])
    array([ 3.,  3.,  3.,  3.,  4.,  5.,  6.,  6.,  6.,  6.])
    
    The same as before, only specify arguments' names for better readability:
    >>> extend_knots_vector(order=4, interval=(3, 6), internal_knots=[4, 5]) 
    array([ 3.,  3.,  3.,  3.,  4.,  5.,  6.,  6.,  6.,  6.])

    >>> extend_knots_vector(4, (3, 6), [4, 5], multiplicities=[3,3])
    array([ 3.,  3.,  3.,  3.,  4.,  4.,  4.,  5.,  5.,  5.,  6.,  6.,  6.,  6.])

    Here we build an extended partition for drawing closed curves:
    >>> extend_knots_vector(order=4, interval=(3, 6), internal_knots=[4, 5], closed=True)
    array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.])

    """
    
    # unpacking the interval extrema
    a, b = interval
    
    # The following could be another way destructure the partition but it is less clear
    # how to distinguish an empty partition (ie. [a, b]) from a filled one.
    #a, *internal_knots, b = internal_knots
    
    if multiplicities is None:
        multiplicities = np.repeat(1, len(internal_knots))

    assert len(internal_knots) == len(multiplicities)

    def repeat_internal_knots_by_multiplicities_in(extended_vector):
        index = order
        for (internal_knot, multiplicity) in zip(internal_knots, multiplicities): 
            extended_vector[index:index+multiplicity] = np.repeat(internal_knot, multiplicity)
            index += multiplicity

    def open_case(): 
        extended_vector = np.empty(order + sum(multiplicities) + order)
        extended_vector[:order] = np.repeat(a, order)
        repeat_internal_knots_by_multiplicities_in(extended_vector)
        extended_vector[-order:] = np.repeat(b, order)

        return extended_vector

    def closed_case(): 
        extended_vector = np.empty(order + sum(multiplicities) + order)
        extended_vector[order-1] = a
        repeat_internal_knots_by_multiplicities_in(extended_vector)
        extended_vector[-order] = b
        
# I'm not able to find a better name to replace `idim`. However
# in the following expression we put evidence on the first two 1s since
# they count knots `a` and `b` respectively.
        idim = 1 + sum(multiplicities) + 1 + (order - 1) 

        for i in range(order-2, -1, -1):
            difference = extended_vector[idim+i-order+1] - extended_vector[idim+i-order]
            extended_vector[i] = extended_vector[i+1] - difference

        for i in range(0, order-1):
            difference = extended_vector[i+order] - extended_vector[i+order-1]
            extended_vector[idim+i] = extended_vector[idim+i-1] + difference

        return extended_vector

    return closed_case() if closed else open_case()
        
def draw(order, interval, internal_knots, control_net, tabs=None, closed=False, multiplicities=None):

    d, n = np.shape(control_net)
    if closed:
        assert n == order + sum(multiplicities) - order + 1
        control_net = np.concatenate((control_net, control_net[:,:order-1]), axis=1)
    else:
        assert n == order + sum(multiplicities)

    if tabs is None:
        a, b = interval
        tabs = np.linspace(start=a, stop=b, num=1000*(b-a))

    extend_knots_vector = extend_knots_partition(
        order, interval, internal_knots, closed, multiplicities)











