# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 17:09:21 2022

@author: phili
"""

import numpy as np
import scipy.special as sp
from InterpolationDuCt import interpolation


from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from BEM_Test import CT_BEM
#plt.close("all")

### Calcul du sillage #########
def fctU(x,y,z,Slicer,SliceValue,Uinf,yawAngle,pitchAngle,rotSpeedR):
    
    #Ct = interpolation(Uinf)
    Ct = CT_BEM(rotSpeedR,pitchAngle,Uinf)[0]
    if (Ct>=0.995):
        Ct = 0.995
    a = 0.5*(1-np.sqrt(1-Ct*(np.cos(yawAngle))**2)) 
    dw = 1+kw*np.log(1+np.exp((x-2*R)/R))
    dv_0 = 0.25*Ct*Uinf*(np.cos(yawAngle))**2*np.sin(yawAngle) #single value
    du_0 = 2*a*Uinf                                      #single value
    dw = 1+kw*np.log(1+np.exp((x-2*R)/R))

    dv = 0.5*dv_0/(dw**2)*(1+sp.erf(x/(R*np.sqrt(2))))
    du = du_0/(dw**2)*0.5*(1+sp.erf(x/(R*np.sqrt(2))))
    
    
    #Calcul de la déflection
    IntSupp = 0   #Valeur de l'intégrale de -infty --> 0 of fctint  = -dv/Uinf 30°
    IntSupp = 0   #Valeur de l'intégrale de -infty --> 0 of fctint  = -dv/Uinf 30°
    
    fctint = -dv/Uinf  #(384)
    fctint2 = fctint[1:]

    yc = np.empty(384)

    dx = 7.89556135770235 # =xv[0][1]  [m]

    #Integral with accumulation of trapezoidal area
    dv[:pD*sl]=0
    yc = dx*(fctint[1:] + fctint[:-1])/2
    yc = np.add.accumulate(yc)

    yc = np.append(yc,yc[-1])
    yc+=IntSupp
    #yc = yc*3
    
    #yc[:pD*sl] = np.nan
    uXYZ = np.zeros((384,96,96))
    
    for i in range(len(x)):
        uXYZ[i] = (Uinf - du[i]*(D**2/(8*sigma_0**2))*np.exp((-(yv3-yc[i]-3*D)**2-(zv3-3*D)**2) /(2*(sigma_0*dw[i])**2)))
        

    if(Slicer == 0):
        # uXYZ[0:int(3*16-8*np.sin(yawAngle)),int(3*16+np.tan(np.pi/2-yawAngle)*(-3*16+hauteur)):96,:] = 9
        return uXYZ[SliceValue,:,:].T, yc
    elif(Slicer== 1):
        return uXYZ[:,SliceValue,:].T, yc
    elif(Slicer== 2):
        return uXYZ[:,:,SliceValue].T, yc
    elif(Slicer == 3):
        return uXYZ[:,:,:], yc
    else:
        print("pas un bon slicer")
        return 0


#fctU(x,y,z,Slicer,hauteur,Uinf,yawAngle):

def fctplot(Slice,SliceValue,yawAngle,pitchAngle,rotSpeedR):
    if(Slice==0):   
        plt.figure(figsize=(6,6))
        uXYZ = fctU(x,y,z,3,SliceValue,Uinf,yawAngle,pitchAngle,rotSpeedR)[0]
        plt.pcolormesh((y)/D,(z)/D,uXYZ[SliceValue,:,:].T/Uinf,cmap='RdBu_r',shading='gouraud');
    
        plt.xlabel('y/D',size = 14)
        plt.ylabel('z/D', size = 14)
        plt.title(r'$\frac{U_{wake}}{U_{inf}}$' r'$\quad$' r'$Y_{aw} =$' +str(yawAngleD), size=20  )
        ax = plt.gca()
        ax.set(xlim=(0,6),ylim=(0,6))

        plt.clim(0.3,1.3);
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        plt.show()
        
    elif(Slice==1):
        plt.figure(figsize=(12,6))
        uXYZ = fctU(x,y,z,3,SliceValue,Uinf,yawAngle,pitchAngle,rotSpeedR)[0]
        plt.pcolormesh((x)/D,(z)/D,uXYZ[:,SliceValue,:].T/Uinf,cmap='RdBu_r',shading='gouraud');
    
        plt.xlabel('x/D',size = 14)
        plt.ylabel('z/D', size = 14)
        plt.title(r'$\frac{U_{wake}}{U_{inf}}$' r'$\quad$' r'$Y_{aw} =$' +str(yawAngleD), size=20  )
        ax = plt.gca()
        ax.set(xlim=(-3,21),ylim=(0,6))
        ax.set_aspect('equal', adjustable='box', anchor='C') #plt.axis('scaled')
        plt.clim(0.2,1.2);
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        plt.show()
        
    elif(Slice == 2):
        plt.figure(figsize=(12,6))
        
        Res = fctU(x,y,z,3,SliceValue,Uinf,yawAngle,pitchAngle,rotSpeedR)
        uXYZ = Res[0]
        yc = Res[1]
        
        plt.pcolormesh(x/D,z/D,uXYZ[:,:,SliceValue].T/Uinf,cmap='RdBu_r',shading='gouraud');
        plt.plot(x/D,yc/D+3,'--')

        plt.xlabel('x/D',size = 14)
        plt.ylabel('y/D', size = 14)
        plt.title(r'$\frac{U_{wake}}{U_{inf}}$' r'$\quad$' r'$Y_{aw} =$' +str(yawAngleD), size=20  )
        ax = plt.gca()
        ax.set(xlim=(-3,21),ylim=(0,6))
        ax.set_aspect('equal', adjustable='box', anchor='C') #plt.axis('scaled')

        plt.clim(0.2,1.2);
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)



D = 126.0
R = D/2
pD = 16   #point par D
nDx = 24   #nD Domaine
nDy = 6
nDz = 6

######### Paramètres simu #########
Uinf = 9 #m/s

kw = 0.044

# Probablement changer cette valeur pour avoir des computes Circle qui reviennenet à Uinf 
#sigma_0 = 0.275*D
sigma_0 = 0.235*D

######### Maillage ######
###
xWT  = 3*D
yWT  = 3*D
zWT  = 3*D

LengthX = nDx*D  #3024 metres
sizeX = nDx*pD  #nX points

LengthY = nDy*D #756 m
sizeY = nDy*pD   #nY points

LengthZ = nDz*D #756 m
sizeZ = nDz*pD   #nZ points


x = np.linspace(0,LengthX,sizeX)
y = np.linspace(0,LengthY,sizeY)
z = np.linspace(0,LengthZ,sizeZ)


x-=xWT
zv3,yv3 = np.meshgrid(y,z)

sl=3
SliceValueX = 5*16
SliceValueY = 3*16
SliceValueZ = 3*16

yawAngleD = 0
yawAngle = yawAngleD*np.pi/180

pitchAngle = 0 
rotSpeed = 10.26  #RPM
rotSpeedR = rotSpeed*2*np.pi/60


Slicer = 0

#fctU(x,y,z,3,SliceValueZ,Uinf,yawAngle,pitchAngle,rotSpeedR)[0]

fctplot(2,SliceValueZ,yawAngle,pitchAngle,rotSpeedR)  #fctplot(Slicer, angle)
