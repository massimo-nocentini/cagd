

import numpy as np

def extend_knots_vector(order, a, b, internal_knots, closed=False, multiplicities=None):
    """
    This function produces an extended knots vector for BSpline functions given a simple one.

    Parameters
    ==========
    a :

    Examples
    ========

    Simple extension for an open curve:
    >>> extend_knots_vector(4, 3, 6, [4, 5])
    array([ 3.,  3.,  3.,  3.,  4.,  5.,  6.,  6.,  6.,  6.])
    
    The same as before, only specify arguments' names for better readability:
    >>> extend_knots_vector(order=4, a=3, b=6, internal_knots=[4, 5]) 
    array([ 3.,  3.,  3.,  3.,  4.,  5.,  6.,  6.,  6.,  6.])

    Here we build an extended partition for drawing closed curves:
    >>> extend_knots_vector(order=4, a=3, b=6, internal_knots=[4, 5], closed=True)
    array([ 0.,  1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.])

    """
    
    # The following could be another way destructure the partition but it is less clear
    # how to distinguish an empty partition (ie. [a, b]) from a filled one.
    #a, *internal_knots, b = internal_knots
    
    if multiplicities is None:
        multiplicities = np.repeat(1, len(internal_knots))

    def repeat_internal_knots_by_multiplicities_in(extended_vector):
        index = order
        for (internal_knot, multiplicity) in zip(internal_knots, multiplicities): 
            extended_vector[index:index+multiplicity] = np.repeat(internal_knot, multiplicity)
            index += multiplicity

    def open_case(): 
        extended_vector = np.empty(order + sum(multiplicities) + order)
        extended_vector[0:order] = np.repeat(a, order)
        repeat_internal_knots_by_multiplicities_in(extended_vector)
        extended_vector[-order:] = np.repeat(b, order)

        return extended_vector

    def closed_case(): 
        extended_vector = np.empty(order + sum(multiplicities) + order)
        extended_vector[order-1] = a
        repeat_internal_knots_by_multiplicities_in(extended_vector)
        extended_vector[-order] = b
        
        idim = 1 + sum(multiplicities) + 1 + (order - 1) 

        for i in range(order-2, -1, -1):
            difference = extended_vector[idim+i-order+1] - extended_vector[idim+i-order]
            extended_vector[i] = extended_vector[i+1] - difference

        for i in range(0, order-1):
            difference = extended_vector[i+order] - extended_vector[i+order-1]
            extended_vector[idim+i] = extended_vector[idim+i-1] + difference

        return extended_vector

    return closed_case() if closed else open_case()
        
