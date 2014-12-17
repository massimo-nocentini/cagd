
import numpy as np
import tensorProductCore as tpc
import tensorProductDrawer as tpd

def exercise_one():
     
    b00 = np.array([0, 0, 0])
    b01 = np.array([2, 0, 0])
    b02 = np.array([4, 0, 0])
    b10 = np.array([0, 2, 0])
    b11 = np.array([2, 2, 0])
    b12 = np.array([4, 2, 2])
    b20 = np.array([0, 4, 0])
    b21 = np.array([2, 4, 4])
    b22 = np.array([4, 4, 4])

    control_net = np.array([[b00,b01,b02],
                            [b10,b11,b12],
                            [b20,b21,b22]], 
                            dtype="float")

    surface = tpc.naive_de_casteljau(control_net)

    tpd.draw(*surface)


#________________________________________________________________________
exercise_one()
