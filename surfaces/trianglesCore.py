
import numpy as np

def u_bar(ntab):
    """
    Produces a set of points uniformly distributed within a given triangle.

    Each point in the produced set is represented in baricentric 
    coordinate form.

    Consumes
    ========
    ntab    number of subinterval in each segment of triangle domain 
            (ie, the parametric domain)

    Produces
    ========
    tri     triangulation matrix of dimension `(ntab^2)x3`, 
            observe that `ntab^2` is the number of triangles within the domain

    U       baricentric coordinates matrix of dimension `3x((ntab+1)*(ntab+2)/2)`
    """

#   observe that it is always possible to halve the next quantity, 
#   since if `ntab` is odd then `ntab+1` is even, and if `ntab` 
#   is even then `ntab+2` is even too, hence both are divisible by 2.
    multi_indeces = (ntab+1)*(ntab+2)/2 

    # points_on_each_segment = ntab+1
    U = np.empty((3, multi_indeces))
    tri = np.empty((ntab**2, 3))

    count = 0
    for kt in range(ntab + 1):
        _np = ntab - kt + 1
        U[:, count:count+_np] = np.array(
            [list(range(_np))[::-1],
             list(range(_np)),
             (kt * np.ones(_np)).tolist()])
        U[:, count:count+_np] /= ntab
        count += _np

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

    return tri, U


def de_casteljau(order, control_net, ntab, V=None):
    """
    Produces a triangular Bezier patch over a baricentric coordinate set.

    For now we don't care about functional case.
    """

    n = order-1
    
    _, Ub = u_bar(n) 

    d, ntot = np.shape(control_net)

    assert ntot == (n+1)*(n+2)/2, "Number of control points mismatch respect to `n`"

    tri, U = u_bar(ntab)
    _, N = np.shape(U)
    surface = np.zeros((d, N, ntot))

    for di in range(d): surface[di,:,:] = np.ones((N,1)) * control_net[di, :ntot]    

    for r in range(n):

        surface_c = np.copy(surface)

        nr = n-r+1

        for k in range(n-r):
            nrk = nr-k
            for i in range(n-r-k):
                ind1 = sum(range(nrk+1, nr+1)) + i
                ind2 = 1 + ind1
                ind3 = sum(range(nrk, nr+1))+ i
                ind  = sum(range(nrk, nr)) + i

                for di in range(d):
                    first_row   = np.multiply(U[0,:], surface_c[di, :, ind1])
                    second_row  = np.multiply(U[1,:], surface_c[di, :, ind2])
                    third_row   = np.multiply(U[2,:], surface_c[di, :, ind3])
                    surface[di,:,ind] = first_row + second_row + third_row

    return surface[:,:,0], tri, U






