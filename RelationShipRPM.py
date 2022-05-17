#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 15:02:05 2022

@author: cheppe
"""
import numpy as np
#from ModeleZongBruno import fctUZongForthGen, ComputeCircle
#WindSpeed = np.linspace(3,25,23)
RotationnalSpeed = np.array([6.972,7.183,7.506,7.942,8.469,9.156,10.926,11.431,11.89,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1,12.1])

def EstimateRotSpeed(Uis):
    if(Uis>12 and Uis<=25):
        Ans =  12.1
    elif(Uis<3):
        print("wind is too slow")
        return np.nan
    elif(Uis>25):
        print("wind is too fast")
        return np.nan
    else:
        x1 = int(Uis)
        x2 = int(Uis+1)
        y1 = RotationnalSpeed[x1-3]
        y2 = RotationnalSpeed[x1-2]
        Ans = y1 +(Uis-x1)*(y2-y1)/(x2-x1)
    return Ans*2*np.pi/60