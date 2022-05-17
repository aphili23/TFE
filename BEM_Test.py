

import os 
os.environ['AFL_PATH'] = "/home/philipin/Documents/source/AirfoilLib"
os.environ['WT4BF_PATH'] = "/home/philipin/Documents/NM_1WT_9ms_LongDom-OPTI_BRU/"
os.environ['BEM_PATH'] = "/home/philipin/Documents/source/uBEM3D"


import uBEM3D 
import numpy as np
from airfoilLib import Airfoil

R = 126.0/2
rho = 1.225
A = np.pi*R**2

def CT_BEM(rotSpeedR,pitchAngle,Uinf):
    af = Airfoil('NREL')
    BEMSet = { 'nBE' : 64 , 'BEtype': '2D' }
    opPnt = {
            'Omega': rotSpeedR,
            'beta' : pitchAngle,
            'V0'   : Uinf}
    BEMSolver = uBEM3D.BEMSolver(af, BEMSet, opPnt)
    T,  M,  Q = BEMSolver.bladeFMCompute(opPnt) #thrust, bending moment, torque
    
    BEM_ret = np.array([T,M,Q])
    CT = T/(0.5*A*rho*Uinf**2)
    return CT, BEM_ret

print("BEM reading...")

#rotSpeed = 10.26   #rpm
#rotSpeedR = rotSpeed*2*np.pi/60
#Uinf = 9.0
#pitchAngle=0.0
# 
#
#print(CT_BEM(rotSpeedR,pitchAngle,Uinf)[1][1]/1000)