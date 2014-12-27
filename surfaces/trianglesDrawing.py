
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def draw(*surfaces, figure_size_tuple=(15,15)):

    sizex, sizey = figure_size_tuple
    matplotlib.rcParams['figure.figsize'] = [sizex, sizey]

    # necessary adjustment if `draw` is used for only one patch
    if len(surfaces) is 2 and not isinstance(surfaces[0], tuple):
        surface, triangles = surfaces
        surfaces = [(surface, triangles)]

    fig = plt.figure() 
    ax = fig.add_subplot(1, 1, 1, projection='3d') 

    for surface, triangles in surfaces:
        x, y, z = surface[0,:],surface[1,:],surface[2,:]
        ax.plot_trisurf(x, y, z, 
                        triangles=triangles, cmap=plt.cm.Spectral)#, edgecolor='none')

    return fig, ax


def repeat_degree_elevation(control_net, snapshots=None, degrees=None, formatting_string="Order {}:")
   
    assert snapshots is not None and degrees is not None, 
        "Or number of snapshots or degrees list must be supplied, but not both."
    
    order, control_net = control_net

    if snapshots:
        def drawer(print_handler):
            runs = 2
            snapshots_list = [int(np.ceil(l)) for l in np.logspace(0,runs,num=snapshots)]
            s = 0
            for i in range(1, (10**runs)+1):
                order, control_net = tc.degree_elevation(order, control_net)
                if i == snapshots[s]:
                    print_handler(order)
                    s += 1

    def print_handler(order):

        if formatting_string is not False: print(formatting_string.format(order))

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        xs = control_net[0,:]
        ys = control_net[1,:]
        zs = control_net[2,:]
        ax.scatter(xs, ys, zs, c='r', marker='o')

        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        plt.show()

    drawer(print_handler) # finally draw some pictures


