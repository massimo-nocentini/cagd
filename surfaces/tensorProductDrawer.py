
from mpl_toolkits.mplot3d import axes3d
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def draw(surface_point_maker, figure_size_tuple=(10,10), squares_per_dim=20, show_figure=True):

    sizex, sizey = figure_size_tuple
    matplotlib.rcParams['figure.figsize'] = [sizex, sizey]

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    squares_per_dim *= 10

    X = np.linspace(0,1,squares_per_dim)
    Y = np.linspace(0,1,squares_per_dim)

    X, Y = np.meshgrid(X, Y)

    surface = np.array([surface_point_maker((u,v)) for u,v in zip(np.ravel(X),np.ravel(Y))])
    sx,sy,sz = surface[:,0],surface[:,1],surface[:,2]

    SX = sx.reshape(X.shape)
    SY = sy.reshape(X.shape)
    SZ = sz.reshape(X.shape)

    surf = ax.plot_surface(SX, SY, SZ, antialiased=True)

    if show_figure: plt.show()

    return surf
