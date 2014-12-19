
import numpy as np
import trianglesCore as tc

def exercise_one():
    control_net = np.array([[0,6,0],
                            [0,3,0],
                            [3,3,6],
                            [0,0,0],
                            [3,0,0], 
                            [6,0,9]]).transpose() 
                            
    surface, tri, U = tc.de_casteljau(3, control_net, 6)

exercise_one()
