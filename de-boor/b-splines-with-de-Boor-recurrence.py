

import numpy as np

def extend_knots_vector(order, a, b, internal_knots, multiplicities=None, closed=False):
    """
    This function produces an extended knots vector for BSpline functions given a simple one.

    Parameters
    ==========
    a:

    Examples
    ========

    Simple extension for an open curve
    >>> extend_knots_vector(4, 3, 6, [4, 5])
    array([ 3.,  3.,  3.,  3.,  4.,  5.,  6.,  6.,  6.,  6.])
    
    The same as before, only specify arguments' names for better readability
    >>> extend_knots_vector(order=4, a=3, b=6, internal_knots=[4, 5]) 
    array([ 3.,  3.,  3.,  3.,  4.,  5.,  6.,  6.,  6.,  6.])

    """
    
    if multiplicities is None:
        multiplicities = np.repeat(1, len(internal_knots))

    def open_case(): 
        extended_vector = np.empty(order + sum(multiplicities) + order)

        extended_vector[0:order] = np.repeat(a, order)

        index = order
        for (internal_knot, multiplicity) in zip(internal_knots, multiplicities): 
            extended_vector[index:index+multiplicity] = np.repeat(internal_knot, multiplicity)
            index += multiplicity
            
        extended_vector[-order:] = np.repeat(b, order)

        return extended_vector

    def closed_case(): ...

    return closed_case() if closed else open_case()
        
