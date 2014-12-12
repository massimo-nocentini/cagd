

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
    
    if multiplicities is None: multiplicities = np.repeat(1, len(internal_knots))

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
        
#       I'm not able to find a better name to replace `idim`. However
#       in the following expression we put evidence on the first two 1s since
#       they count knots `a` and `b` respectively.
        idim = 1 + sum(multiplicities) + 1 + (order - 1) 

        for i in range(order-2, -1, -1):
            difference = extended_vector[idim+i-order+1] - extended_vector[idim+i-order]
            extended_vector[i] = extended_vector[i+1] - difference

        for i in range(0, order-1):
            difference = extended_vector[i+order] - extended_vector[i+order-1]
            extended_vector[idim+i] = extended_vector[idim+i-1] + difference

        return extended_vector

    return closed_case() if closed else open_case()
        

def de_Boor(extended_knots_partition, control_net, tabs): 
    """
    Produces BSpline curve interpolating the given control net.
    """
    
    if type(tabs) is not np.ndarray: tabs = np.array(tabs)

    control_net = control_net.transpose()

    d, n = np.shape(control_net)
    order = len(extended_knots_partition) - n

    C = np.zeros((d, len(tabs)))

    def foreach_dimension(f):
        for i in range(d): f(i)

    ind = 0 

    for r in range(order-1, n):
        
        sup_extrema = extended_knots_partition[r+1]
        tloc = tabs[np.where(extended_knots_partition[r] <= tabs)]
        if r is n-1:
            tloc = tloc[np.where(tloc <= sup_extrema)]
        else:
            tloc = tloc[np.where(tloc < sup_extrema)]

        nloc = len(tloc)

#       The following handle the case where no tabulation point belong to 
#       the corrent subinterval between consecutive internal knots.
#       This could happen if a non uniform tabulation set of points is given.
        if not nloc: continue

        Qloc = np.zeros((order, nloc, d))

#       The following function can be seen as the very first de Boor's algorithm step:
#       for each dimension and for each tabulation point, it `copies` the relative
#       control point component, the same for each tabulation point. Therefore 
#       we can see it as the first step, ie. the initializer one.
        def build_temp_matrix(i):
            spline_support = np.matrix(control_net[i,r-order+1:r+1]).transpose()
            Qloc[:,:,i] = spline_support.dot(np.ones((1, nloc)))
        foreach_dimension(build_temp_matrix)

#       `j` starts from 1 since the `0` iteration is done by the initialization above.
#       In total, the outer `for` runs `order-1` times.
        for j in range(1, order):
            alfa = np.zeros((order-j, nloc)) 
#           For a fixed `j`, the inner `for` runs `order-j` times, therefore the
#           comprehensive complexity is (order-1)*order/2, in other words O(order^2)
            for i in range(0, order-j):
#               The following shows because for a change of very first or very last knots 
#               (or both) there isn't a change in curve's shape: they never get in `alfa`.
                inf_extrema = extended_knots_partition[i+1+r-order+j]
                sup_extrema = extended_knots_partition[i+1+r]
                knots_distance = sup_extrema - inf_extrema
                alfa[i,:] = (tloc - inf_extrema) / knots_distance if knots_distance > 0 else 0

                def baricentric_combine(s): 
                    Qloc[i,:,s] = (1-alfa[i,:])*Qloc[i,:,s] + alfa[i,:]*Qloc[i+1,:,s] 
                foreach_dimension(baricentric_combine)

        def fill_in_curve_points(i): C[i, ind:ind+nloc] = Qloc[0, :, i]
        foreach_dimension(fill_in_curve_points)

        ind += nloc

    assert np.shape(C) == np.shape(C[:, :ind])

    return C.transpose()


def knot_insertion(t_hat, extended_knots_partition, control_net, order):
    """
    Produces an augmented knots partition and control net.
    """

    n, _ = np.shape(control_net)
    augmented_knots_partition = np.zeros(1+len(extended_knots_partition))
    
    r = None 
#   Before the inf extrema of interval, namely `a`, it's forbidden
#   It is also forbidden beyond the sup extrema, namely `b`.
    for r in range(order-1, n+1): 
        if t_hat < extended_knots_partition[r]: break

    augmented_knots_partition[:r] = extended_knots_partition[:r] 
    augmented_knots_partition[r] = t_hat
    augmented_knots_partition[r+1:] = extended_knots_partition[r:] 

    r -= 1 # index of inf extrema of interval containing `t_hat`

    def omega(i, s=order):
        knots_slack = extended_knots_partition[i+s-1]-extended_knots_partition[i] 
        hat_difference = t_hat - extended_knots_partition[i]
        return hat_difference / knots_slack if knots_slack > 0 else 0

    combination_matrix = np.zeros((n+1,n))

    for i in range(0, r-order+2):   combination_matrix[i,i] = 1
    for i in range(r-order+2, r+1): o = omega(i); combination_matrix[i,i-1:i+1] = [1-o, o] 
    for i in range(r+1, n+1):   combination_matrix[i,i-1] = 1

    augmented_control_net = combination_matrix.dot(control_net)

    return augmented_knots_partition, augmented_control_net


def raise_internal_knots_to_max_smooth(
        order, internal_knots, multiplicities, extended_vector, 
        control_net, fix_extended_vector, closed=False):
    """
    Generator that raise each internal knot to max smoothness.

    It keeps applying the knot insertion algorithm on each internal knot 
    until its multiplicity gets `order-2`, the max request for smoothness.

    It is a generator of tuples, where each tuple consists of multiplicities, 
    extended_vector and control_net at each rising step. Observe that 
    internal knots partition argument isn't returned since it doesn't 
    change because only multiplicities of knots already there are raised.
    """

    if closed: control_net = np.concatenate((control_net, control_net[:order-1, :]), axis=0)

    for knot, m in zip(range(len(internal_knots)), range(len(multiplicities))):

        while multiplicities[m] < order-2:

            extended_vector, control_net = knot_insertion(
                internal_knots[knot], extended_vector, control_net, order)

            multiplicities[m] += 1

            if closed:
#               First fix the extended partion to restore cyclic property
                extended_vector = fix_extended_vector(multiplicities)

#               Second, fix control net tail part to be a cyclic repetition of the head
                if knot != len(internal_knots)-1:
                    control_net[-(order-1):, :] = control_net[:order-1, :]
                else:
                    control_net[:order-1, :] = control_net[-(order-1):, :]

            yield multiplicities, extended_vector, control_net

