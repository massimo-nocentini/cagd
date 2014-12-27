
import numpy as np

def u_bar(ntab, return_multi_indices_matrix=False, triangles_partitions=False):
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

    In the following we show how to partition triangles:
    >>> tri, U, partitioned_triangles = u_bar(ntab=5, triangles_partitions=True)
    >>> partitioned_triangles['upside']
    [(0, 1, 6), (1, 2, 7), (2, 3, 8), (3, 4, 9), (4, 5, 10), (6, 7, 11), (7, 8, 12), (8, 9, 13), (9, 10, 14), (11, 12, 15), (12, 13, 16), (13, 14, 17), (15, 16, 18), (16, 17, 19), (18, 19, 20)]

    >>> partitioned_triangles['upside_down']
    [(1, 6, 7), (2, 7, 8), (3, 8, 9), (4, 9, 10), (7, 11, 12), (8, 12, 13), (9, 13, 14), (12, 15, 16), (13, 16, 17), (16, 18, 19)]

    >>> partitioned_triangles['on_left_inv_diagonal']
    [(0, 1, 6), (1, 2, 7), (2, 3, 8), (3, 4, 9), (4, 5, 10)]

    >>> partitioned_triangles['on_right_diagonal']
    [(0, 1, 6), (6, 7, 11), (11, 12, 15), (15, 16, 18), (18, 19, 20)]

    >>> partitioned_triangles['on_bottom_diagonal']
    [(4, 5, 10), (9, 10, 14), (13, 14, 17), (16, 17, 19), (18, 19, 20)]

    """

#   observe that it is always possible to halve the next quantity, 
#   since if `ntab` is odd then `ntab+1` is even, and if `ntab` 
#   is even then `ntab+2` is even too, hence both are divisible by 2.
    multi_indeces = int((ntab+1)*(ntab+2)/2) 

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

    # the following lists allow to partition triangles
    partitioned_triangles = {
        'upside':[],
        'upside_down':[],
        'on_left_inv_diagonal':[],
        'on_right_diagonal':[],
        'on_bottom_diagonal':[]
        }

    def update_tri_matrix(a, b, c):
        update_tri_matrix.count += 1
        tri[update_tri_matrix.count,:] = np.array([a, b, c])

    update_tri_matrix.count = -1

    for kt in range(ntab-1):

        nk = ntab+2-kt
        sm = sum(range(nk,ntab+2))
        end = sm + (ntab-kt-1)

        for i, ind in enumerate(range(sm, end)):

            upside_triangle = (ind, ind+1, ind+nk-1)
            upside_down_triangle = (ind+1, ind+nk-1, ind+nk)

            update_tri_matrix(*upside_triangle)
            update_tri_matrix(*upside_down_triangle)
            
            partitioned_triangles['upside'].append(upside_triangle) 
            partitioned_triangles['upside_down'].append(upside_down_triangle) 

#           using `i` from the enumeration allow us to look for the very first
#           triangle without comparing against `sm`, the start value of `range`
            if i is 0: partitioned_triangles['on_right_diagonal'].append(upside_triangle) 

        last_triangle = (end, end+1, end+nk-1)
        update_tri_matrix(*last_triangle)
        partitioned_triangles['upside'].append(last_triangle) 
        partitioned_triangles['on_bottom_diagonal'].append(last_triangle) 

    rightmost_bottom_triangle = (multi_indeces-3, multi_indeces-2, multi_indeces-1)
    update_tri_matrix(*rightmost_bottom_triangle)
    partitioned_triangles['upside'].append(rightmost_bottom_triangle) 
    partitioned_triangles['on_right_diagonal'].append(rightmost_bottom_triangle) 
    partitioned_triangles['on_bottom_diagonal'].append(rightmost_bottom_triangle) 

    partitioned_triangles['on_left_inv_diagonal'] = partitioned_triangles['upside'][:ntab]

    assert update_tri_matrix.count == (ntab**2 - 1)

    assert (len(partitioned_triangles['on_left_inv_diagonal']) ==
            len(partitioned_triangles['on_right_diagonal']) ==
            len(partitioned_triangles['on_bottom_diagonal']) == 
            ntab)

    result = (tri, U)
    if return_multi_indices_matrix: result += (multi_indices_matrix,)
    if triangles_partitions: result += (partitioned_triangles,)
    
    return result

#_______________________________________________________________________

def de_casteljau(order, control_net, ntab=None, triangulation=None, subdivision=False):
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

    if ntab is None: ntab = 2*n # double the number of control points for "tabulation"

    if triangulation is None: 
        tri, U, multi_indices_matrix = u_bar(ntab, return_multi_indices_matrix=True)
    else:
#       we unpack with a collecting `rest`...
        tri, U, *rest = triangulation 
        
#       ...and do some dimensional checks
        _n, N = np.shape(U)
        T, t = np.shape(tri)
        assert N == (ntab+1)*(ntab+2)/2 and _n is 3
        assert T == ntab**2 and t is 3

    _, N = np.shape(U)
    surface = np.zeros((d, N, ntot))

    def initialize_surface_with_control_points(di): 
        surface[di,:,:] = np.ones((N,1)) * control_net[di, :ntot]    

    foreach_dimension(initialize_surface_with_control_points)
   
    if subdivision: 

#       Prepare subpatches: everyone has the same shape of the original control net
        left_subpatch = np.empty(np.shape(control_net))
        right_subpatch = np.empty(np.shape(control_net))
        bottom_subpatch = np.empty(np.shape(control_net))
        
#       Indices vectors (aka. `diagonals`) for sub-patches construction
#       Pay attention: it is important to not use `ntab` here, since
#       `ntab` deals with tabulation only, ie only with the set of triangles
#       we'll use to represent patches, which can extend arbitrarily.
        offsets = np.array(range(order)) 
        left_inv_diagonal = np.array(range(order))
        right_diagonal = np.cumsum([0] + list(range(order,1,-1)))
        bottom_diagonal = np.cumsum([n] + list(range(n,0,-1)))

#       Vector-wise initialization with control points on main diagonals
        left_subpatch[:, left_inv_diagonal] = control_net[:, left_inv_diagonal]
        right_subpatch[:, right_diagonal] = control_net[:, right_diagonal]
        bottom_subpatch[:, bottom_diagonal] = control_net[:, bottom_diagonal]

    layers = []

    for r in range(n):

        surface_c = np.copy(surface)
        nr = n-r+1

        layer = []
        layers.append(layer)

        for k in range(n-r):

            nrk = nr-k

            for i in range(n-r-k):

#               Indices of each triangle extrema to combine in this step
                ind1 = sum(range(nrk+1, nr+1)) + i
                ind2 = 1 + ind1
                ind3 = sum(range(nrk, nr+1)) + i

#               Index of `new` position for combined point
                ind  = sum(range(nrk, nr)) + i

#               Combined points using indices and barycentric coordinates
#               in each triangle needed for the combination.
#               Observe that for indexing in the barycentri coordinates dimension
#               we refer to `ntab`, hence we're independent from the choice
#               of barycentric coordinates (aka. "tabulation")
                top_point = (surface[:, 0, ind1] 
                            + surface[:, ntab, ind1] 
                            + surface[:, -1, ind1])/float(3)

                left_point = (  surface[:, 0, ind2] 
                                + surface[:, ntab, ind2] 
                                + surface[:, -1, ind2])/float(3)

                right_point = ( surface[:, 0, ind3] 
                                + surface[:, ntab, ind3] 
                                + surface[:, -1, ind3])/float(3)

#               de Casteljau interpolation point for the current iteration
#               Observe that `if r == n-1` then this point __really lies on patch__.
                combined_point = (top_point + left_point + right_point)/float(3)

                layer.append(combined_point)

                def update_surface(di):
                    first_row   = np.multiply(U[0,:], surface_c[di, :, ind1])
                    second_row  = np.multiply(U[1,:], surface_c[di, :, ind2])
                    third_row   = np.multiply(U[2,:], surface_c[di, :, ind3])
                    surface[di,:,ind] = first_row + second_row + third_row

                foreach_dimension(update_surface)

#   At the end of layers construction and if subdivision is required we build
#   control nets relative to left, right and bottom sub-patches respectively.
    if subdivision:       
        
        for r, layer in zip(range(n), map(np.array, layers)):

#           Remember: """AttributeError: 'numpy.ndarray' object has no attribute 'pop'"""

            left_inv_diagonal = left_inv_diagonal[1:]
            left_inv_diagonal += offsets[-(r+1)]
            left_subpatch[:, left_inv_diagonal] = layer[
                offsets[:n-r]].transpose()

            right_diagonal = right_diagonal[1:]
            right_diagonal -= offsets[-1:-order+r:-1]
            right_subpatch[:, right_diagonal] = layer[
                np.cumsum([0] + list(range(n-r,1,-1)))].transpose()

            bottom_diagonal = bottom_diagonal[1:] 
            bottom_diagonal -= offsets[-1:-order+r:-1] + 1
            bottom_subpatch[:, bottom_diagonal] = layer[
                np.cumsum([n-1-r] + list(range(n-1-r,0,-1)))].transpose()

    surface = surface[:,:,0]

#   Making results tuple, not matter the optional requests about triangulation
#   and subdivision, we return a nested tuple in both cases.
    results = ((surface, layers),)
    if triangulation is None: results += ((tri, U, multi_indices_matrix),)
    if subdivision: results += ((left_subpatch, right_subpatch, bottom_subpatch),)

    return results

#_______________________________________________________________________

def degree_elevation(order, control_net):
    """
    Produces an augmented control net to interpolate the same surface but with increased degree.
    """

    d, ntot = np.shape(control_net)
    old_order = order
    new_order = old_order + 1

    _, U = u_bar(ntab=old_order)

    tri, _, partitioned_triangles = u_bar(ntab=old_order-1, triangles_partitions=True)

    _, N = np.shape(U)

    sandbox = np.zeros((d, N, ntot))

    augmented_control_net = np.empty((d, new_order))

    def foreach_dimension(func):
        for di in range(d): func(di)

    def initialize_sandbox_with_control_points(di): 
        sandbox[di,:,:] = np.ones((N,1)) * control_net[di,:]    

    foreach_dimension(initialize_sandbox_with_control_points)
  
    covered_points = []

    def rotate_clockwise(triangle): 
        a, b, c = triangle
        return a, c, b

    def log(*args):
        log = """
control point position {} is under assignment 
from combination coordinate {} 
using upside down triangle {} (rotated clockwise {})
        """
        print(log.format(*args))

    def upside_down_triangles_handler(triangles):
        t = 0
        for top_most_vertex_in_right_diagonal_triangles, forward_offset in zip(
                np.cumsum([old_order + 2] + list(range(old_order,old_order-2 -1,-1))),
                range(old_order-2, 0, -1)): 

            print("start of diagonal iteration on upside down triangles for control point position {}, with forward offset {}:".
                    format(top_most_vertex_in_right_diagonal_triangles, forward_offset))

            for k, comb in zip([top_most_vertex_in_right_diagonal_triangles + i 
                            for i in range(forward_offset)],
                         np.cumsum([old_order+1+forward_offset] + list(range(old_order-1,old_order-forward_offset,-1)))): 

                log(k, comb, triangles[t], rotate_clockwise(triangles[t]))

                covered_points.append(k)
                t += 1
          
    def on_left_diagonals_triangles_handler(triangles):           
        print("on LEFT diagonal:")            
        for l, triangle in zip(range(1, old_order), triangles):
            comb = old_order - l 
            log(l, comb, triangle, rotate_clockwise(triangle))
            covered_points.append(l)

    def on_right_diagonals_triangles_handler(triangles):
        print("on RIGHT diagonal:")            
        right_diagonal = [r for r in np.cumsum([old_order+1] + list(range(old_order,2,-1)))]
        for (ri, r), triangle in zip(enumerate(right_diagonal), triangles):
            comb = right_diagonal[-(ri+1)] 
            log(r, comb, triangle, rotate_clockwise(triangle))
            covered_points.append(r)

    def on_bottom_diagonals_triangles_handler(triangles):
        print("on BOTTOM diagonal:")            
        bottom_diagonal = [b for b in np.cumsum([2*old_order] + list(range(old_order-1,1,-1)))]
        for (bi, b), triangle in zip(enumerate(bottom_diagonal), triangles):
            comb = bottom_diagonal[-(bi+1)]
            log(b, comb, triangle, rotate_clockwise(triangle))
            covered_points.append(b)

    upside_down_triangles_handler(partitioned_triangles['upside_down']) 
    on_left_diagonals_triangles_handler(partitioned_triangles['on_left_inv_diagonal'])
    on_right_diagonals_triangles_handler(partitioned_triangles['on_right_diagonal'])
    on_bottom_diagonals_triangles_handler(partitioned_triangles['on_bottom_diagonal'])

    covered_points.append(0)
    covered_points.append(old_order)
    covered_points.append(int((old_order+1)*(old_order+2)/2)-1)
    
    assert sorted(covered_points) == list(range(int((old_order+1)*(old_order+2)/2)))
    print("assertion satisfied: all new points are covered")

    return new_order, augmented_control_net



