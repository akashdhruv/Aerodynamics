import numpy as np
import matplotlib.pyplot as plt
from math import *

# Defining Grid
N = 200 
xStart,xEnd = -4.0,4.0
yStart,yEnd = -2.0,2.0 
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)    
X,Y = np.meshgrid(x,y)          

# Defining Freestream
Uinf = 1.0
alphaInDegrees = 0.0
alpha = alphaInDegrees*pi/180

# Calculating Freestream components
uFreestream = Uinf*cos(alpha)
vFreestream = Uinf*sin(alpha)

# Calculating Freestream Stream-Function
psiFreestream = + Uinf*cos(alpha)*Y - Uinf*sin(alpha)*X

# Defining function to calculate vortex properties
def flowvortex(ga,xv,yv,X,Y):
    u = +ga/(2*pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = -ga/(2*pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    psi = ga/(4*pi)*np.log((X-xv)**2+(Y-yv)**2)
    return u,v,psi

# Defining row of vortices    
Nv=500
xv=np.linspace(-100,100,Nv)
yv=np.zeros(np.size(xv))
ga=.001

# Defining velocity vectors for the Grid
uvortex=np.zeros((N,N),dtype=float)
vvortex=np.zeros((N,N),dtype=float)
psivortex=np.zeros((N,N),dtype=float)

# Calculating Induced Velocities
for i in range(Nv):
    u1,v1,psi1,=flowvortex(ga,xv[i],yv[i],X,Y)
    uvortex=uvortex+u1
    vvortex=vvortex+v1
    psivortex=psivortex+psi1
    

# Plotting streamlines
size = 10
plt.figure()
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uvortex,vvortex,density=4.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xv,yv,c='#CD2305',s=80,marker='o')
plt.contour(X,Y,psivortex,levels=[0.0],colors='#CD2305',linewidths=2,linestyles='--')
plt.show()  

