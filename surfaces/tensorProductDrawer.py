
from mpl_toolkits.mplot3d import axes3d
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def draw(surface_point_maker, figure_size_tuple=(10,10), squares_per_dim=20, show_figure=True):

    sizex, sizey = figure_size_tuple
    matplotlib.rcParams['figure.figsize'] = [sizex, sizey]


    squares_per_dim *= 10

    X = np.linspace(0,1,squares_per_dim)
    Y = np.linspace(0,1,squares_per_dim)

    X, Y = np.meshgrid(X, Y)

    surface = np.array([surface_point_maker((u,v)) for u,v in zip(np.ravel(X),np.ravel(Y))])
    sx,sy,sz = surface[:,0],surface[:,1],surface[:,2]

    SX = sx.reshape(X.shape)
    SY = sy.reshape(X.shape)
    SZ = sz.reshape(X.shape)

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
