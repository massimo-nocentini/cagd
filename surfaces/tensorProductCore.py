
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

def naive_de_casteljau(control_net, tabs=None, squares_per_dim=20):

    squares_per_dim *= 10

    X = np.linspace(0,1,squares_per_dim)
    Y = np.linspace(0,1,squares_per_dim)

    X, Y = np.meshgrid(X, Y)

    surface = np.array([de_casteljau_surface((u,v), control_net) 
                            for u,v in zip(np.ravel(X),np.ravel(Y))])

    return surface, X.shape

def vectorized_de_casteljau(control_net, tabs=None, breaks=100):
    """
    Vectorized implementation of de Casteljau algorithm.
    """

    if tabs is None: 
        tabs = np.linspace(0,1, breaks), np.linspace(0,1, breaks)

    u_tab, v_tab = tabs

#   Sestini's code calls `m` what we call `rows` and
#   calls `n` what we call `cols`.
    rows, cols, d = np.shape(control_net)

    k = min(rows, cols)
    mtab, ntab = len(u_tab), len(v_tab)
    u_tab, v_tab = np.meshgrid(u_tab, v_tab) 
    shape = u_tab.shape
    print(u_tab, u_tab.shape)
    #u_tab, v_tab = u_tab.reshape((mtab, 1)), v_tab.reshape((ntab, 1))
    u_tm, v_tm = 1-u_tab, 1-v_tab

    surface = np.zeros((mtab, ntab, d, rows, cols))

    def initialize(it, jt, i_d): surface[it, jt, i_d, :, :] = control_net[:, :, i_d]
    [initialize(it, jt, i_d) for it in range(mtab) for jt in range(ntab) for i_d in range(d)]

    def update(i_d, i, j): 
        uv_tm_tm = np.multiply(u_tm, v_tm)
        uv_tm_tab = np.multiply(u_tm, v_tab)
        uv_tab_tm = np.multiply(u_tab, v_tm)
        uv_tab_tab = np.multiply(u_tab, v_tab)
        surface[:, :, i_d, i, j] = (
                np.multiply(uv_tm_tm, surface[:, :, i_d, i, j])
                + np.multiply(uv_tm_tab, surface[:, :, i_d, i, j+1])
                + np.multiply(uv_tab_tm, surface[:, :, i_d, i+1, j])
                + np.multiply(uv_tab_tab, surface[:, :, i_d, i+1, j+1]))
    
    [update(i_d, i, j)  for i_d in range(d)
                        for r in range(1, k) 
                            for i in range(rows-r) 
                                for j in range(cols-r) ]

    return surface[:, :, :, 0, 0], shape
