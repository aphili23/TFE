# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 11:10:42 2021

@author: portable
"""



import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

from RelationShipRPM import EstimateRotSpeed 

from ModeleBruno import fctU, Uinf, kw
from BEM_Test import CT_BEM

# WT specific packages
import readField as rf
from readDiskData import readDiskData

import os



D = 126.0
R = D/2

pD = 16   #point par D
nDx = 24   #nD Domaine
nDy = 6
nDz = 6

Uw = 0


sl = 3 #Placement de la 1ere eolienne


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
LengthX = nDx*D  #3024 metres
sizeX = nDx*pD  #nX points

def ComputeCircle(Uyz_CC,yawAngleI): #Renvoie la moyenne sur 1 cercle #yawAngleI pour la WT i où on compute circle
                                    
    count = 0
    tot = 0
    for i in range(len(z)):
        for j in range(len(y)):
            if(np.cos(yawAngleI)<=1*10**(-5)):
                return 0
            elif(((z[i]-3*D)/R)**2+((y[j]-3*D)/(R*np.cos(yawAngleI)))**2)<=1.0:
                count+=1
                tot += Uyz_CC[i,j]

    return tot/count



def fctUZongForthGen(NbrEo,dist,yawAngleTab,pitchAngleTab,rotSpeedTabR,Uinf,ControlRotSpeed):
    Somme = 0
    Ud = 0
    Uis = 0
    ycTab = []
    UisTab = np.zeros(3)
    if(ControlRotSpeed[0]==False):
        rotSpeedR = EstimateRotSpeed(Uinf) #[rad/s]
        rotSpeedTabR[0] = rotSpeedR
    for i in range(NbrEo):
        if(i==0):      #Première Eo #sl -1 : car on fait computecircle 1 D avant l'éolienne.
            Uis = Uinf
            Somme += (Uis-fctU(x-dist*D*i,y,z,3,0,Uis,yawAngleTab[i],pitchAngleTab[i],rotSpeedTabR[i])[0])**2 #slice 3 because 3d dimention
            Ud = ComputeCircle(Uinf - np.sqrt(Somme[(sl-1+dist*(i+1))*16,:,:].T),yawAngleTab[i+1]) #Ud pour la 1ere eo
            yc = fctU(x-dist*D*i,y,z,3,0,Uis,yawAngleTab[i],pitchAngleTab[i],rotSpeedTabR[i])[1]
            yc[:pD*(sl+i*dist)] = np.nan
            ycTab.append(yc)
            UisTab[i] = Uinf 
            UisTab[i+1] = Ud
            if(ControlRotSpeed[i+1]==False):
                rotSpeedR = EstimateRotSpeed(Ud) #[rad/s]
                rotSpeedTabR[i+1] = rotSpeedR
                
        elif(i==NbrEo-1): #Dernière Eo
            Uis = Ud
            Somme += (Uis-fctU(x-dist*D*i,y,z,3,0,Uis,yawAngleTab[i],pitchAngleTab[i],rotSpeedTabR[i])[0])**2 #slice 3 because 3d dimention
            yc = fctU(x-dist*D*i,y,z,3,0,Uis,yawAngleTab[i],pitchAngleTab[i],rotSpeedTabR[i])[1]
            yc[:pD*(sl+i*dist)] = np.nan
            ycTab.append(yc)
            
        else: 
            Uis = Ud
            Somme += (Uis-fctU(x-dist*D*i,y,z,3,0,Uis,yawAngleTab[i],pitchAngleTab[i],rotSpeedTabR[i])[0])**2 #slice 3 because 3d dimention
            Ud = ComputeCircle(Uinf - np.sqrt(Somme[(sl-1+dist*(i+1))*16,:,:].T),yawAngleTab[i+1]) #Ud pour la 1ere eo
            yc = fctU(x-dist*D*i,y,z,3,0,Uis,yawAngleTab[i],pitchAngleTab[i],rotSpeedTabR[i])[1]
            yc[:pD*(sl+i*dist)] = np.nan
            ycTab.append(yc)
            UisTab[i+1] = Ud
            if(ControlRotSpeed[i+1]==False):
                rotSpeedR = EstimateRotSpeed(Ud) #[rad/s]
                rotSpeedTabR[i+1] = rotSpeedR
                
    Uw = Uinf - np.sqrt(Somme)
    return Uw,ycTab,UisTab

def fctplot(Slicer,SliceValue,yawAngleTab,pithAngleTab,rotSpeedTabR,Uinf,NbrEo,dist):
    if(Slicer==0):
        plt.figure(figsize=(6,6))
        uXYZ,ycTab,UisTab = fctUZongForthGen(NbrEo,dist,yawAngleTab,pitchAngleTab,rotSpeedTabR,Uinf,ControlRotSpeed)
      
        plt.pcolormesh((y)/D,(z)/D,uXYZ[SliceValue,:,:].T/Uinf,cmap='RdBu_r',shading='gouraud');
                
        # plt.plot(x/D,fctU(x,y,z,2,90,Uinf)[1]/D,'--')
        plt.xlabel('y/D',size = 14)
        plt.ylabel('z/D', size = 14)
        #plt.title(r'$\frac{U_{wake}}{U_{inf}}$' r'$\quad$' r'$Y_{aw} =$' +str(yawAngle), size=20)
        ax = plt.gca()
        ax.set(xlim=(0,6),ylim=(0,6))
       
        # colorbar
        plt.clim(0.3,1.3);
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
    elif(Slicer==1):
        
        plt.figure(figsize=(12,6))
        
        uXYZ,ycTab,UisTab = fctUZongForthGen(NbrEo,dist,yawAngleTab,pitchAngleTab,rotSpeedTabR,Uinf,ControlRotSpeed)
        plt.pcolormesh((x)/D,(z)/D,uXYZ[:,SliceValue,:].T/Uinf,cmap='RdBu_r',shading='gouraud');
    
        plt.xlabel('x/D',size = 14)
        plt.ylabel('z/D', size = 14)
        #plt.title(r'$\frac{U_{wake}}{U_{inf}}$' r'$\quad$' r'$Y_{aw} =$' +str(yawAngle), size=20)
        ax = plt.gca()
        ax.set(xlim=(-3,21),ylim=(0,6))
        ax.set_aspect('equal', adjustable='box', anchor='C') #plt.axis('scaled')

       
        # colorbar
        plt.clim(0.2,1.2);
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)

    elif(Slicer==2):
                
        plt.figure(figsize=(12,6))
        uXYZ,ycTab,UisTab = fctUZongForthGen(NbrEo,dist,yawAngleTab,pitchAngleTab,rotSpeedTabR,Uinf,ControlRotSpeed)

        plt.pcolormesh((x)/D,(y)/D,uXYZ[:,:,SliceValue].T/Uinf,cmap='RdBu_r',shading='gouraud');
        for i in range(NbrEo):
            plt.plot((x-dist*i)/D,ycTab[i]/D+3,'--')
        # plt.plot(x/D,fctU(x,y,z,2,90,Uinf)[1]/D,'--')
        plt.xlabel('x/D',size = 14)
        plt.ylabel('y/D', size = 14)
        #plt.title(r'$\frac{U_{wake}}{U_{inf}}$' r'$\quad$' r'$Y_{aw} =$' +str(yawAngle), size=20)
        ax = plt.gca()
        ax.set(xlim=(-3,21),ylim=(0,6))
        ax.set_aspect('equal', adjustable='box', anchor='C') #plt.axis('scaled')

        # colorbar
        plt.clim(0.2,1.2);
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(cax=cax)
        

def ComputePuissance(uXYZ):   #Calcul de la puissance grace à la BEM
    puissanceTab = np.zeros(NbrEo)
    UisTab = np.zeros(NbrEo)
    for i in range(NbrEo):
        sl1 = uXYZ[2*16+i*dist*16,:,:].T
        Uis = ComputeCircle(sl1,yawAngleTab[i])
        torque = CT_BEM(rotSpeedTabR[i],pitchAngleTabR[i],Uis)[1][2]
        puissanceTab[i] = torque*rotSpeedTabR[i]
        UisTab[i] = Uis
    return puissanceTab/10**6,UisTab   #[MW]~~


######### DEBUT PARAMETRES DE SIMU #########

SliceValueX = 4*16
SliceValueY = 3*16
SliceValueZ = 3*16


dist = 7
NbrEo = 3

yawAngleTab = np.array([0.0,0.0,0.0])    
yawAngleTabR = yawAngleTab*np.pi/180

pitchAngleTab = np.array([0.0,0.0,0.0])    #deg
pitchAngleTabR = pitchAngleTab*np.pi/180   #rad

rotSpeedTab = np.array([7.8889,8.15,8.15])   #rpm
rotSpeedTab = np.array([10.26,10.26,10.26])   #rpm
rotSpeedTabR = rotSpeedTab*2*np.pi/60       #rad/s

RangeRotSpeed  = np.linspace(5,12.1,10)   #rpm
RangeRotSpeedR = RangeRotSpeed*2*np.pi/60  #rad/s
RangeRotSpeedR_rev = RangeRotSpeedR[::-1]

PuissanceTab_WT0 = np.empty(len(RangeRotSpeed))
PuissanceTab_WT1 = np.empty(len(RangeRotSpeed))
PuissanceTab_WT2 = np.empty(len(RangeRotSpeed))
ControlRotSpeed = [False,False,False] #False -> L'éolienne calcule elle meme la rot speed

PuissanceTab_WT_all = np.empty(len(RangeRotSpeed))



#uXYZ = fctUZongForthGen(NbrEo,dist,yawAngleTabR,pitchAngleTabR,rotSpeedTabR,Uinf,ControlRotSpeed)[0]

######### FIN PARAMETRES DE SIMU #########



######### DEBUT - Plot de 2 slices différentes #########
"""
fctplot(2,SliceValueZ,yawAngleTab,pitchAngleTab,rotSpeedTabR,Uinf,NbrEo,dist)
fctplot(0,SliceValueX,yawAngleTab,pitchAngleTab,rotSpeedTabR,Uinf,NbrEo,dist)
"""
######### FIN - Plot de 2 slices différentes #########

   

######## CTRL - DEBUT OPTI RANGESPEED 1WT & 2 WT #########
"""
#Verifier que le control est activé
ControlRotSpeed = [True,True,False] #False -> L'éolienne calcule elle meme la rot speed

PuissanceTab_WT0 = np.empty(((len(RangeRotSpeed)),len(RangeRotSpeed)))
PuissanceTab_WT1 = np.empty(((len(RangeRotSpeed)),len(RangeRotSpeed)))
PuissanceTab_WT2 = np.empty(((len(RangeRotSpeed)),len(RangeRotSpeed)))

PuissanceTab_WT_all = np.empty(((len(RangeRotSpeed)),len(RangeRotSpeed)))
UisTabXY = []

X,Y = np.meshgrid(RangeRotSpeed,RangeRotSpeed)
for i in range(len(RangeRotSpeedR)):
    for j in range(len(RangeRotSpeedR)):
        rotSpeedTabR[0] = RangeRotSpeedR[i]
        rotSpeedTabR[1] = RangeRotSpeedR_rev[j]
        uXYZ_pw = fctUZongForthGen(NbrEo,dist,yawAngleTabR,pitchAngleTabR,rotSpeedTabR,Uinf,ControlRotSpeed)[0]
        Puiss = ComputePuissance(uXYZ_pw)
        UisTabXY.append(Puiss[1]) #Renvoie un tableau de taille 3 (Uis 1,2,3)
        PuissTab = Puiss[0]
        PuissanceTab_WT0[i][j]  = PuissTab[0]
        PuissanceTab_WT1[i][j]  = PuissTab[1]
        PuissanceTab_WT2[i][j]  = PuissTab[2]
        PuissanceTab_WT_all[i][j] = sum(PuissTab)
        print("Etat : "+ str(len(RangeRotSpeedR)*i+j))


MaxPuissFarm = np.max(PuissanceTab_WT_all)
Index_MaxPuissFarm = np.where(PuissanceTab_WT_all==np.max(PuissanceTab_WT_all))   #Indice dans la valeur P_max ds le tableau PuissMax
UisTab_MAX_Power = UisTabXY[Index_MaxPuissFarm[0][0]*len(RangeRotSpeed)+Index_MaxPuissFarm[1][0]]

#f= open("Result.txt","w+")
#f.truncate()
#f.write("Results OPTI WT0 and WT1 for U_inf = " + str(Uinf) + "\n \n")
#f.write("The max power reach is " + str(round(MaxPuissFarm,3)) + " [MW]" + "\n")
#f.write("Velocities Uis  are " + str(UisTab_MAX_Power)+ "[m/s] \n")
#f.write("The rot speed WT0_OPT is " + str(RangeRotSpeed[Index_MaxPuissFarm[0][0]]) + " RPM \n")
#f.write("The rot speed WT1_OPT is " + str(RangeRotSpeed[Index_MaxPuissFarm[1][0]]) + " RPM \n")
#f.close()

#PLot 3D de la puissance, avec controle de rotSpeed des WT1 & WT 2.
fig = plt.figure(1,figsize=(10,10))
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, PuissanceTab_WT_all, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.title("Power of the farm influenced by the rotationspeed \n of the " + "$1^{st}$" + " and " "$2^{nd}$" + " wind turbine",size=20)
ax.set_xlabel(r'$\omega_{WT_1} \: [rpm] $',size=18)
ax.set_ylabel(r'$\omega_{WT_2} \: [rpm] $',size=18)
ax.set_zlabel('Power [MW]',size=20)
#plt.savefig("figures/Opti3D_Omega1&2WT.png")

#PLot 2D de la puissance, avec controle de rotSpeed de WT1. RotSpeed de WT 2 optimal. 
plt.figure(2,figsize=(7.5,5.5))
plt.plot(RangeRotSpeed,PuissanceTab_WT_all[:,Index_MaxPuissFarm[-1][0]]) #coupe selon y
plt.title("Slice of the 3D Power \n " + '$\omega_{WT_2,OPT}$ = ' + str(round(RangeRotSpeed[Index_MaxPuissFarm[-1][0]],2)) + " RPM",size=18)
plt.xlabel(r'$\omega_{WT_1} \: [rpm]$',size=20)
plt.ylabel('Power [MW]',size=20)
plt.savefig("figures/Opti3D_sliceWT2_fixed.png")
plt.show()


#PLot 2D de la puissance, avec controle de rotSpeed de WT2. RotSpeed de WT1 optimal. 
plt.figure(3,figsize=(7.5,5.5))
plt.plot(RangeRotSpeed,PuissanceTab_WT_all[Index_MaxPuissFarm[0][0],:]) #coupe selon y
plt.title("Slice of the 3D Power \n " + '$\omega_{WT_1,OPT}$ = ' + str(round(RangeRotSpeed[Index_MaxPuissFarm[0][0]],2)) + " RPM",size=18)
plt.xlabel(r'$\omega_{WT_2} \: [rpm]$',size=20)
plt.ylabel('Power [MW]',size=20)
#plt.savefig("figures/Opti3D_sliceWT1_fixed.png")

plt.show()
"""
######## CTRL - FIN OPTI RANGESPEED 1WT & 2 WT #########



######## CTRL - DEBUT OPTI RANGESPEED 1WT #########
"""
#Verifier que le control 1WT est activé
ControlRotSpeed = [True,False,False] #False -> L'éolienne calcule elle meme la rot speed
for i in range(len(RangeRotSpeedR)):
    rotSpeedTabR[0] = RangeRotSpeedR[i]
    uXYZ_pw = fctUZongForthGen(NbrEo,dist,yawAngleTabR,pitchAngleTabR,rotSpeedTabR,Uinf,ControlRotSpeed)[0]
    PuissTab = ComputePuissance(uXYZ_pw)[0]
    PuissanceTab_WT0[i]  = PuissTab[0]
    PuissanceTab_WT1[i]  = PuissTab[1]
    PuissanceTab_WT2[i]  = PuissTab[2]
    PuissanceTab_WT_all[i] = sum(PuissTab)
    
plt.figure(figsize=(9,5))
plt.plot(RangeRotSpeed,PuissanceTab_WT0,color='red')
plt.plot(RangeRotSpeed,PuissanceTab_WT1,color='green')
plt.plot(RangeRotSpeed,PuissanceTab_WT2,color='blue')
plt.plot(RangeRotSpeed,PuissanceTab_WT_all,color='black')

plt.title("Power of the 3 WT influenced by \n the rotation speed of the 1st WT",size=18)

plt.xlabel(r'$\omega_{WT_1} \: [rpm]$',size=20)
plt.ylabel('Power [MW]',size=20)
ax = plt.gca()
ax.legend(['WT0','WT1','WT2','SUM'])

#Placement d'un point avec coord
PuissanceTab_WT0 = np.around(PuissanceTab_WT0,3)
RangeRotSpeed = np.around(RangeRotSpeed,3)

P_MAX_WT0 = max(PuissanceTab_WT0)
index = np.where(PuissanceTab_WT0==P_MAX_WT0)
Omega_MAX_WT0 = RangeRotSpeed[index[0]]

plt.text(float(Omega_MAX_WT0),float(P_MAX_WT0-0.3),(round(Omega_MAX_WT0[0],1),round(P_MAX_WT0,2)))
plt.plot(Omega_MAX_WT0[0],P_MAX_WT0,'r*')

#Placement d'un point avec coord
PuissanceTab_WT1 = np.around(PuissanceTab_WT1,3)

P_MAX_WT1 = max(PuissanceTab_WT1)
index = np.where(PuissanceTab_WT1==P_MAX_WT1)
Omega_MAX_WT1 = RangeRotSpeed[index[0]]

plt.text(Omega_MAX_WT1-0.2,P_MAX_WT1-0.4,(round(Omega_MAX_WT1[0],1),round(P_MAX_WT1,2)))
plt.plot(Omega_MAX_WT1[0],P_MAX_WT1,'g*')


#Placement d'un point avec coord
PuissanceTab_WT2 = np.around(PuissanceTab_WT2,3)

P_MAX_WT2 = max(PuissanceTab_WT2)
index = np.where(PuissanceTab_WT2==P_MAX_WT2)
Omega_MAX_WT2 = RangeRotSpeed[index[0]]

plt.text(Omega_MAX_WT2+0.1,P_MAX_WT2-0.4,(round(Omega_MAX_WT2[0],1),round(P_MAX_WT2,2)))
plt.plot(Omega_MAX_WT2[0],P_MAX_WT2,'b*')

#Placement d'un point avec coord
PuissanceTab_WT_all = np.around(PuissanceTab_WT_all,3)

P_MAX_WT_all = max(PuissanceTab_WT_all)
index = np.where(PuissanceTab_WT_all==P_MAX_WT_all)
Omega_MAX_WT_all = RangeRotSpeed[index[0]]

plt.text(Omega_MAX_WT_all+0.1,P_MAX_WT_all-0.4,(round(Omega_MAX_WT_all[0],1),round(P_MAX_WT_all,2)))
plt.plot(Omega_MAX_WT_all[0],P_MAX_WT_all,'k*')


plt.savefig("figures/Opti_Omega1WT.png")
"""
######## CTRL - FIN OPTI RANGESPEED 1WT #########





######## CTRL - DEBUT OPTI RANGESPEED 2WT #########
"""
#Verifier que le control 2WT est activé
ControlRotSpeed = [False,True,False] #False -> L'éolienne calcule elle meme la rot speed
for i in range(len(RangeRotSpeedR)):
    print("Etat : " +str(i))
    rotSpeedTabR[1] = RangeRotSpeedR[i]
    uXYZ_pw = fctUZongForthGen(NbrEo,dist,yawAngleTabR,pitchAngleTabR,rotSpeedTabR,Uinf,ControlRotSpeed)[0]
    PuissTab = ComputePuissance(uXYZ_pw)[0]
    PuissanceTab_WT0[i]  = PuissTab[0]
    PuissanceTab_WT1[i]  = PuissTab[1]
    PuissanceTab_WT2[i]  = PuissTab[2]
    PuissanceTab_WT_all[i] = sum(PuissTab)
    
plt.figure(figsize=(9,5))
plt.plot(RangeRotSpeed,PuissanceTab_WT0,color='red')
plt.plot(RangeRotSpeed,PuissanceTab_WT1,color='green')
plt.plot(RangeRotSpeed,PuissanceTab_WT2,color='blue')
plt.plot(RangeRotSpeed,PuissanceTab_WT_all,color='black')

plt.title("Power of the 3 WT influenced by \n the rotation speed of the 2nd WT",size=18)

plt.xlabel(r'$\omega_{WT_2} \: [rpm]$',size=20)
plt.ylabel('Power [MW]',size=20)
ax = plt.gca()
ax.legend(['WT0','WT1','WT2','SUM'])

##Placement d'un point avec coord
#PuissanceTab_WT0 = np.around(PuissanceTab_WT0,3)
#RangeRotSpeed = np.around(RangeRotSpeed,3)
#
#P_MAX_WT0 = max(PuissanceTab_WT0)
#index = np.where(PuissanceTab_WT0==P_MAX_WT0)
#Omega_MAX_WT0 = RangeRotSpeed[index[0]]
#
#plt.text(float(Omega_MAX_WT0),float(P_MAX_WT0-0.3),(round(Omega_MAX_WT0[0],1),round(P_MAX_WT0,2)))
#plt.plot(Omega_MAX_WT0[0],P_MAX_WT0,'r*')

#Placement d'un point avec coord
#PuissanceTab_WT1 = np.around(PuissanceTab_WT1,6)
#
#P_MAX_WT1 = max(PuissanceTab_WT1)
#index = np.where(PuissanceTab_WT1==P_MAX_WT1)
#Omega_MAX_WT1 = RangeRotSpeed[index[0]]
#
#plt.text(Omega_MAX_WT1-0.2,P_MAX_WT1-0.4,(round(Omega_MAX_WT1[0],1),round(P_MAX_WT1,2)))
#plt.plot(Omega_MAX_WT1[0],P_MAX_WT1,'g*')


#Placement d'un point avec coord
#PuissanceTab_WT2 = np.around(PuissanceTab_WT2,3)
#
#P_MAX_WT2 = max(PuissanceTab_WT2)
#index = np.where(PuissanceTab_WT2==P_MAX_WT2)
#Omega_MAX_WT2 = RangeRotSpeed[index[0]]
#
#plt.text(Omega_MAX_WT2+0.1,P_MAX_WT2-0.4,(round(Omega_MAX_WT2[0],1),round(P_MAX_WT2,2)))
#plt.plot(Omega_MAX_WT2[0],P_MAX_WT2,'b*')
#
##Placement d'un point avec coord
#PuissanceTab_WT_all = np.around(PuissanceTab_WT_all,3)
#
#P_MAX_WT_all = max(PuissanceTab_WT_all)
#index = np.where(PuissanceTab_WT_all==P_MAX_WT_all)
#Omega_MAX_WT_all = RangeRotSpeed[index[0]]
#
#plt.text(Omega_MAX_WT_all+0.1,P_MAX_WT_all-0.4,(round(Omega_MAX_WT_all[0],1),round(P_MAX_WT_all,2)))
#plt.plot(Omega_MAX_WT_all[0],P_MAX_WT_all,'k*')


#plt.savefig("figures/Opti_Omega2WT.png")
"""
######## CTRL - FIN OPTI RANGESPEED 2sWT #########




############# DEBUT DE CALIBRATION #############

#Verifier que ctrl = false :> Calibration - R
ControlRotSpeed = [False,False,False] #False -> L'éolienne calcule elle meme la rot speed
#ControlRotSpeed = [True,True,True] #False -> L'éolienne calcule elle meme la rot speed

rotSpeedTab = np.array([10.296,10.296,10.296])   #rpm
rotSpeedTabR = rotSpeedTab*2*np.pi/60       #rad/s

uXYZ = fctUZongForthGen(NbrEo,dist,yawAngleTabR,pitchAngleTabR,rotSpeedTabR,Uinf,ControlRotSpeed)[0]

#Paramètres de simu BF

simuName = "NREL_9ms_TI06"
# for fields:
D    = 126.0;
Uinf = 9.0
xWT  = 3*D
yWT  = 3*D
zWT  = 3*D
# for disk data
nt  = 12000



#Décalage du x pour etre identique à celui de Zong.
x+=xWT

## Plot slice BF à x = xSliceValue
##############################################
#X,U = rf.readStats('%s_avg'%(simuName));
#y,z,Us1 = rf.getSlice(X, U, 0, x[SliceValueX]); #slicing is necessary for stats, 7D après
#
#plt.figure(figsize=(8,8))
#plt.pcolormesh((y-yWT)/D,z/D,Us1/Uinf,cmap='RdBu_r',shading='gouraud');
#plt.xlabel('y/D',size =20)
#plt.ylabel('z/D',size =20)
##plt.title(r'$\frac{U_{wake}}{U_{inf}}$' r'$\quad -$' r'$\quad$' r'$Yaw =$' +str(angle) + '°' + r'$\quad -$' r'$\quad$' r'$U_{inf} = \: 9 \: m/s$', size=20)
#ax = plt.gca()
#ax.set(xlim=(-5,5),ylim=(-2,8))
#        
## draw WT line
#plt.plot([0,0],[0,6],'--',color="blue")
#plt.plot([-3,3],[3,3],'--',color="blue")
#
## colorbar
#plt.clim(0.3,1.2);
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#cbar = plt.colorbar(cax=cax)
#cbar.set_label("      " + r'$\frac{U_{wake}}{U_{inf}}$',rotation=360,size=35)
#plt.show()
##############################################



X,U = rf.readStats('%s_avg'%(simuName));

MeanErrorVel = np.empty(len(x))
MeanVelBF = np.empty(len(x))
MeanVelZong = np.empty(len(x))

MeanErrorVel[-1]=np.nan
MeanErrorVel[-2]=np.nan

MeanVelBF[-1]=np.nan
MeanVelBF[-2]=np.nan

MeanVelZong[-1]=np.nan
MeanVelZong[-2]=np.nan


for i in range(len(x)-2):
    y,z,Us = rf.getSlice(X, U, 0,x[i]);
    U_meanBF = ComputeCircle(Us,0)
    U_meanZong = ComputeCircle(uXYZ[i,:,:],0)

    MeanVelBF[i] = U_meanBF
    MeanVelZong[i] = U_meanZong
    MeanErrorVel[i] = U_meanBF-U_meanZong

#Graphe de Mean Error Velocity in % along x. Difference avec computeCirlce des 2 modèles
    
plt.figure(figsize=(16,8))
MeanErrorVel[-1] = np.nan
Error = MeanErrorVel/Uinf*100
plt.plot(x/D,Error)
plt.grid("on")
#plt.xlabel(r'$\frac{x}{D}$',size=20)
#plt.ylabel(r'$\frac{U_{BF}-U_{Zong}}{U_{\infty}}$',size=20)
#plt.title("Mean error in percentage of \n the velocity along x",size=24)
#plt.plot([-3,21],[0,0],'--')
#plt.plot([-1,-1],[-10,10],'--','r')
#plt.plot([6,6],[-10,10],'--','r')
#plt.plot([13,13],[-10,10],'--','r')

#ax = plt.gca()
#ax.set(xlim=(-3,21),ylim=(-15,40))
#ax.set_aspect('equal', adjustable='box', anchor='C') #plt.axis('scaled')

#Graphe de Mean velocity BF - computeCirle along X for BF
plt.savefig("Mean error MODIF for k_wall = " + str(kw) + ".png")
plt.figure(figsize=(16,8))
plt.plot(x/D,MeanVelBF)

#Graphe de Mean Velocity Zong - computeCirle along X for Zong
plt.figure(5,figsize=(16,8))
plt.plot(x/D,MeanVelZong)

#Vérification de l'erreur pour différentes valeurs de x (commence à 0)
ind_ToCheck = 16*np.array([2,2+dist,2+2*dist])

#
for i in range(3):
    print("##################")
    print("u_is Zong vaut " +str(ComputeCircle(uXYZ[ind_ToCheck[i],:,:],0)))
    y,z,Us = rf.getSlice(X, U, 0,x[ind_ToCheck[i]]);
    print("Le x de BF vaut : " +str(x[ind_ToCheck[i]]/D))
    print("u_is BF vaut " +str(ComputeCircle(Us,0)))
    print("Error in percent at x = " +str(ind_ToCheck[i]/16) +" D is : " +str(Error[ind_ToCheck[i]]))

############# FIN DE CALIBRATION #############


