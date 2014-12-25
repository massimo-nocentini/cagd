
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
