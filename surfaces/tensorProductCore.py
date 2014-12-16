
import numpy as np

def de_casteljau(t, control_net):
    """
    Returns a point on the Bezier curve.

    This function applies de Casteljau algorithm, the
    constructive approach (ie, the geometric one).

    Given a control net, it produces a tuple (p, ud, ld) such that:
        - p is the point on the Bezier curve that interpolate the control net;
        - ud is an array containing points on the upper diagonal for curve splitting;
        - ld is an array containing points on the lower diagonal for curve splitting.

    The overall complexity is O(n^2), where n is the number of control points.
    """

    n, d = dimensions = np.shape(control_net)

    Q = np.copy(control_net)

    upper_diagonal = np.empty(dimensions)
    lower_diagonal = np.empty(dimensions)

    upper_diagonal[0,:] = Q[0,:]
    lower_diagonal[0,:] = Q[-1,:]
    
    for k in range(1,n):
        for i in range(n-k):
            Q[i,:] = (1-t)*Q[i,:] + t*Q[i+1,:]
        
        upper_diagonal[k,:] = Q[0,:]
        lower_diagonal[k,:] = Q[-(k+1),:]
    
    return Q[0,:], upper_diagonal, lower_diagonal


def de_casteljau_surface(params, control_net):
    """
    Produces a point on the Bezier surface.

    Examples
    ========

    The following example works on a square control net:
    >>> b00 = np.array([0, 0, 0])
    >>> b01 = np.array([2, 0, 0])
    >>> b02 = np.array([4, 0, 0])
    >>> b10 = np.array([0, 2, 0])
    >>> b11 = np.array([2, 2, 0])
    >>> b12 = np.array([4, 2, 2])
    >>> b20 = np.array([0, 4, 0])
    >>> b21 = np.array([2, 4, 4])
    >>> b22 = np.array([4, 4, 4])
    >>> control_net = np.array([[b00,b01,b02],\
                                [b10,b11,b12],\
                                [b20,b21,b22]], dtype="float")
    >>> u,v = .5, .5
    >>> point = de_casteljau_surface((u,v), control_net)
    >>> point
    array([ 2.,  2.,  1.])
    """

    control_net = np.copy(control_net)

    u, v = params

    rows, cols, d = np.shape(control_net)
    
    def combine_col(r, c): return (1-v)*control_net[r, c, :] + v*control_net[r, c+1, :]
    def combine_row(r, c): return (1-u)*control_net[r, c, :] + u*control_net[r+1, c, :]

    # for now assume to work with a square control net
    for k in range(rows-1):
        for r in range(rows-k): 
            for c in range(cols-k-1): control_net[r, c, :] = combine_col(r, c) 

        # after the previous loops we've eliminated the rightmost column
        for c in range(cols-k-1):
            for r in range(rows-k-1): control_net[r, c, :] = combine_row(r, c)

    return control_net[0, 0, :]
