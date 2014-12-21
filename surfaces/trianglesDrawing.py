
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def draw(surface, triangles, figure_size_tuple=(15,15)):

    sizex, sizey = figure_size_tuple
    matplotlib.rcParams['figure.figsize'] = [sizex, sizey]

    x, y, z = surface[0,:],surface[1,:],surface[2,:]

    fig = plt.figure() 
    ax = fig.add_subplot(1, 1, 1, projection='3d') 
    ax.plot_trisurf(x, y, z, triangles=triangles, cmap=plt.cm.Spectral)#, edgecolor='none')
