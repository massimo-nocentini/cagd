
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

    multi_indeces = (ntab+1)*(ntab+2)/2
    points_on_each_segment = ntab+1

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
