
from mpl_toolkits.mplot3d import axes3d
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def draw(surface, shape, figure_size_tuple=(15,15), show_figure=True):

    sizex, sizey = figure_size_tuple

    matplotlib.rcParams['figure.figsize'] = [sizex, sizey]
    
    sx,sy,sz = surface[:,0],surface[:,1],surface[:,2]
    SX = sx.reshape(shape)
    SY = sy.reshape(shape)
    SZ = sz.reshape(shape)

    surf = drawer(lambda ax: ax.plot_surface(SX, SY, SZ, antialiased=True), show_figure)

#   Maybe the following could be factored better about `zdir`
    xcontour = drawer(lambda ax: ax.contour(SX, SY, SZ, zdir='x', antialiased=True), show_figure)
    ycontour = drawer(lambda ax: ax.contour(SX, SY, SZ, zdir='y', antialiased=True), show_figure)
    zcontour = drawer(lambda ax: ax.contour(SX, SY, SZ, zdir='z', antialiased=True), show_figure)

    return surf, xcontour, ycontour, zcontour

def drawer(plotter, show_figure=True):
    """
    Produce a plot on a newly created figure, showing it if desired.
    """
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    plot = plotter(ax)
    if show_figure: plt.show()
    return plot
