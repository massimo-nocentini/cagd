
import numpy as np
import tensorProductCore as tpc
import tensorProductDrawer as tpd

def exercise_two():
     
    b00 = np.array([0, 0, 0])
    b01 = np.array([2, 0, 0])
    b02 = np.array([4, 0, 0])
    b03 = np.array([6, 0, 0])
    b10 = np.array([0, 2, 0])
    b11 = np.array([2, 2, 4])
    b12 = np.array([4, 2, 4])
    b13 = np.array([6, 2, 0])
    b20 = np.array([0, 4, 0])
    b21 = np.array([2, 4, 4])
    b22 = np.array([4, 4, 4])
    b23 = np.array([6, 4, 0])
    b30 = np.array([0, 6, 0])
    b31 = np.array([2, 6, 0])
    b32 = np.array([4, 6, 0])
    b33 = np.array([6, 6, 0])

    control_net = np.array([[b00,b01,b02, b03],
                            [b10,b11,b12, b13],
                            [b20,b21,b22, b23], 
                            [b30,b31,b32, b33]], 
                            dtype="float")

    surface_maker = lambda params: tpc.de_casteljau_surface(params, control_net)

    tpd.draw(surface_maker)


#________________________________________________________________________
exercise_two()
