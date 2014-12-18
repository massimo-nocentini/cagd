
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
    tri = np.empty((ntab^2, 3))

    count = 0
    for kt in range(ntab + 1):
        np = ntab - kt + 1
        U[:, count:count+np] = np.array(
            [list(range(np-1,1,-1)),
             list(range(np)),
             kt * np.ones(np)], dtype="f")
        U[:, count:count+np] /= ntab
        count += np

    count = -1
    for kt in range(ntab-1):

        nk = ntab+2-kt
        sm = sum(range(nk,ntab+2))
        ind = None

        for it in range(ntab-kt-1):
            ind = sm+it+1 # this `+1` is suspect
            count += 1
            tri[count,:] = np.array([ind, ind+1, ind+nk-1])
            count += 1
            tri[count,:] = np.array([ind+1, ind+nk-1, ind+nk])

        count += 1
        tri[count,:] = np.array([ind+1, ind+2, ind+nk])

    count += 1
#   Maybe in the following indeces we should subtract 1 if they are used to index matrix U
    tri[count,:] = np.array([multi_indeces-2, multi_indeces-1, multi_indeces])

    assert count == ntab^2 - 1

    return tri, U


def de_casteljau(order, control_net, ntab, V=None):
    """
    Produces a triangular Bezier patch over a baricentric coordinate set.

    For now we don't care about functional case.
    """

    n = order-1
    
    trib, Ub = u_bar(n) 

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
                ind1 = sum(range(nrk+1, nr+1)) + i + 1
                ind2 = 1 + ind1
                ind3 = sum(range(nrk, nr+1))+ i + 1
                ind  = sum(nrk, nr) + i + 1

                for di in range(d):
                    first_row   = np.multiply(U[1,:], surface_c[di, :, ind1])
                    second_row  = np.multiply(U[2,:], surface_c[di, :, ind2])
                    third_row   = np.multiply(U[3,:], surface_c[di, :, ind3])
                    surface[di,:,ind] = first_row + second_row + third_row

    return surface[:,:,0], tri, U






