
import numpy as np

def u_bar(ntab, return_multi_indices_matrix=False):
    """
    Produces a set of points uniformly distributed within a given triangle.

    Each point in the produced set is represented in baricentric 
    coordinate form. Also, a set of triangles is produced and, if desired,
    the matrix of multi indices is produced too.

    Observe that triangles set is quite different from the one produced
    by Matplotlib (which does use Delaunay algorithm?), and it is based 
    on Prof. Sestini algorithm.

    Consumes
    ========
    ntab    number of subinterval in each side of triangle domain 
            (ie, the parametric domain). Observe that on each tringle line
            will lie `ntab+1` points.

    return_multi_indices_matrix     
            `True` will include multi indices matrix in tuple result.

    Produces
    ========
    tri     triangulation matrix of dimension `(ntab^2)x3`, 
            observe that `ntab^2` is the number of triangles within the domain.

    U       baricentric coordinates matrix of dimension `3x((ntab+1)*(ntab+2)/2)`
            observe that `((ntab+1)*(ntab+2)/2)` is the number of different
            multi indices required to make a base.

    multi_indices_matrix
            multi indices matrix of internal points (ie, triangles extrema), 
            returned if argument `return_multi_indices_matrix` is `True`.

    Examples
    ========
    >>> tri, U, multi_indices = u_bar(ntab=4, return_multi_indices_matrix=True)
    >>> multi_indices
    array([[ 4.,  3.,  2.,  1.,  0.,  3.,  2.,  1.,  0.,  2.,  1.,  0.,  1.,
             0.,  0.],
           [ 0.,  1.,  2.,  3.,  4.,  0.,  1.,  2.,  3.,  0.,  1.,  2.,  0.,
             1.,  0.],
           [ 0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,  1.,  2.,  2.,  2.,  3.,
             3.,  4.]])
    >>> U
    array([[ 1.  ,  0.75,  0.5 ,  0.25,  0.  ,  0.75,  0.5 ,  0.25,  0.  ,
             0.5 ,  0.25,  0.  ,  0.25,  0.  ,  0.  ],
           [ 0.  ,  0.25,  0.5 ,  0.75,  1.  ,  0.  ,  0.25,  0.5 ,  0.75,
             0.  ,  0.25,  0.5 ,  0.  ,  0.25,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.25,  0.25,  0.25,  0.25,
             0.5 ,  0.5 ,  0.5 ,  0.75,  0.75,  1.  ]])
    >>> tri
    array([[  0.,   1.,   5.],
           [  1.,   5.,   6.],
           [  1.,   2.,   6.],
           [  2.,   6.,   7.],
           [  2.,   3.,   7.],
           [  3.,   7.,   8.],
           [  3.,   4.,   8.],
           [  5.,   6.,   9.],
           [  6.,   9.,  10.],
           [  6.,   7.,  10.],
           [  7.,  10.,  11.],
           [  7.,   8.,  11.],
           [  9.,  10.,  12.],
           [ 10.,  12.,  13.],
           [ 10.,  11.,  13.],
           [ 12.,  13.,  14.]])
    """

#   observe that it is always possible to halve the next quantity, 
#   since if `ntab` is odd then `ntab+1` is even, and if `ntab` 
#   is even then `ntab+2` is even too, hence both are divisible by 2.
    multi_indeces = (ntab+1)*(ntab+2)/2 

    U = np.empty((3, multi_indeces))
    tri = np.empty((ntab**2, 3))

    count = 0
    for kt in range(ntab + 1):
        _np = ntab - kt + 1
        U[:, count:count+_np] = np.array(
            [list(range(_np))[::-1],
             list(range(_np)),
             (kt * np.ones(_np)).tolist()])
        count += _np

    multi_indices_matrix = np.copy(U) # just have a copy of multi indices
    U /= ntab # make the matrix represent baricentric coordinates

    def update_tri_matrix(a, b, c):
        update_tri_matrix.count += 1
        tri[update_tri_matrix.count,:] = np.array([a, b, c])

    update_tri_matrix.count = -1

    for kt in range(ntab-1):

        nk = ntab+2-kt
        sm = sum(range(nk,ntab+2))
        end = sm + (ntab-kt-1)

        for ind in range(sm, end):
            update_tri_matrix(ind, ind+1, ind+nk-1)
            update_tri_matrix(ind+1, ind+nk-1, ind+nk)

        update_tri_matrix(end, end+1, end+nk-1)

    update_tri_matrix(multi_indeces-3, multi_indeces-2, multi_indeces-1)

    assert update_tri_matrix.count == (ntab**2 - 1)

    return (tri, U, multi_indices_matrix) if return_multi_indices_matrix else (tri, U)


def de_casteljau(order, control_net, ntab, triangulation=None):
    """
    Produces a triangular Bezier patch over a baricentric coordinate set.

    For now we don't care about functional case, for more about that a porting
    of Sestini's implementation is pending.

    Consumes
    ========
    order   Order of the surface

    control_net     
            Matrix of points of the form `3x((n+1)*(n+2)/2)` where `n = order-1`

    ntab    Number of subinterval in each side of the domain triangle

    triangulation
            an iterable such that the first element is a set of triangles,
            and the second element is a matrix of baricentric coordinates, after
            the iterable can contain anything, which is ignored.
            If `None`, a triangulation pair is obtained using `u_bar` function,
            in this case that computed pair is included in the returned result tuple.    

    Produces
    ========
    surface 
            A matrix of points on the triangular Bezier patch, with form 
            `d x ((ntab+1)*(ntab+2)/2)`, the same as baricentric coordinate matrix
    """

    n = order-1
    
    d, ntot = np.shape(control_net)

    assert ntot == (n+1)*(n+2)/2, "Number of control points mismatch respect to `n`"

    def foreach_dimension(func):
        for di in range(d): func(di)

    if triangulation is None: 
        tri, U, multi_indices = u_bar(ntab, return_multi_indices_matrix=True)
    else:
        tri, U, *rest = triangulation # we unpack with a collecting `rest`

    _, N = np.shape(U)
    surface = np.zeros((d, N, ntot))

    def initialize_surface_with_control_points(di): 
        surface[di,:,:] = np.ones((N,1)) * control_net[di, :ntot]    

    foreach_dimension(initialize_surface_with_control_points)

    for r in range(n):

        surface_c = np.copy(surface)
        nr = n-r+1

        for k in range(n-r):

            nrk = nr-k

            for i in range(n-r-k):

                ind1 = sum(range(nrk+1, nr+1)) + i
                ind2 = 1 + ind1
                ind3 = sum(range(nrk, nr+1)) + i
                ind  = sum(range(nrk, nr)) + i

                def update_surface(di):
                    first_row   = np.multiply(U[0,:], surface_c[di, :, ind1])
                    second_row  = np.multiply(U[1,:], surface_c[di, :, ind2])
                    third_row   = np.multiply(U[2,:], surface_c[di, :, ind3])
                    surface[di,:,ind] = first_row + second_row + third_row

                foreach_dimension(update_surface)

    surface = surface[:,:,0]

    return (surface, tri, U, multi_indices_matrix) if triangulation is None else surface
