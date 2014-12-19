
import numpy as np
import trianglesCore as tc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def exercise_one():
   
    control_net_from_Farin = [  [0,6,0],
                                [0,3,0],
                                [3,3,6], 
                                [0,0,0], 
                                [3,0,0], 
                                [6,0,9]]  

    control_net = np.array(control_net_from_Farin[::-1]).transpose() 

    surface, tri, U = tc.de_casteljau(3,control_net, 10) 
    
    x, y, z = surface[0,:],surface[1,:],surface[2,:]
    
    figsize(15,15) 
    fig = plt.figure() 
    ax = fig.add_subplot(1, 1, 1, projection='3d') 
    ax.plot_trisurf(x, y, z, cmap=plt.cm.Spectral)#, edgecolor='none')

#________________________________________________________________________
exercise_one()
