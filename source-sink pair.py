import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 200                           # Number of points in each direction
xStart,xEnd = -4.0,4.0            # x-direction boundaries
yStart,yEnd = -2.0,2.0            # y-direction boundaries
x = np.linspace(xStart,xEnd,N)    # x 1D-array
y = np.linspace(yStart,yEnd,N)    # y 1D-array
X,Y = np.meshgrid(x,y)            # generation of the mesh grid 

Uinf = 1.0        # freestream speed
alphaInDegrees = 0.0       # angle of attack (in degrees)
alpha = alphaInDegrees*pi/180

# computing the velocity components on the mesh grid
uFreestream = Uinf*cos(alpha)
vFreestream = Uinf*sin(alpha)

# computing the stream-function on the mesh grid
psiFreestream = + Uinf*cos(alpha)*Y - Uinf*sin(alpha)*X

def flowproperties(strength,xs,ys,X,Y):
    u = strength/(2*pi)*(X-xs)/((X-xs)**2+(Y-ys)**2)
    v = strength/(2*pi)*(Y-ys)/((X-xs)**2+(Y-ys)**2)
    psi = strength/(2*pi)*np.arctan2((Y-ys),(X-xs))
    return u,v,psi
def flowdoublet(mu,xs,ys,X,Y):
    u = - mu/(2*pi)*((X-xd)**2-(Y-yd)**2)/((X-xd)**2+(Y-yd)**2)**2
    v = - mu/(2*pi)*2*(X-xd)*(Y-yd)/((X-xd)**2+(Y-yd)**2)**2
    psi = -mu/(2*pi)*(Y-yd)/((X-xd)**2+(Y-yd)**2)
    return u,v,psi

def flowvortex(ga,xv,yv,X,Y):
    u = +ga/(2*pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = -ga/(2*pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    psi = ga/(4*pi)*np.log((X-xv)**2+(Y-yv)**2)
    return u,v,psi

""" ssource=5.0

xsource,ysource=-1.0,0.0

usource,vsource,psisource=flowproperties(ssource,xsource,ysource,X,Y)

ssink=-5.0
xsink,ysink=1.0,0.0
usink,vsink,psisink=flowproperties(ssink,xsink,ysink,X,Y)

u = uFreestream + usource + usink
v = vFreestream + vsource + vsink
psi = psiFreestream + psisource +psisink """

xd,yd=0.0,0.0
mu=5
ga=5
xv,yv,=0.0,0.0
udoublet,vdoublet,psidoublet=flowdoublet(mu,xd,yd,X,Y)
uvortex,vvortex,psivortex=flowvortex(ga,xv,yv,X,Y)

u=udoublet+uFreestream+uvortex
v=vdoublet+vFreestream+vvortex
psi=psidoublet+psiFreestream+psivortex

# plotting

size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xd,yd,c='#CD2305',s=80,marker='o')

""" # computing the stagnation point
xStagnation = xsource - ssource/(2*pi*Uinf)*cos(alpha)
yStagnation = ysource - ssource/(2*pi*Uinf)*sin(alpha)

xstag=xsink - ssink/(2*pi*Uinf)*cos(alpha)
ystag=ysink - ssink/(2*pi*Uinf)*sin(alpha)

# adding the stagnation point to the figure
plt.scatter([xStagnation,xstag],[yStagnation,ystag],c='b',s=80,marker='o') """

# adding the dividing line to the figure
#if (alpha==0.0):
plt.contour(X,Y,psi,linewidths=2,linestyles='--')

plt.show()

Cp = 1.0-(u**2+v**2)/Uinf**2

plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
"""plt.scatter([xsource,xsink],[ysource,ysink],c='#CD2305',s=80,marker='o')

plt.scatter([xStagnation,xstag],[yStagnation,ystag],c='b',s=80,marker='o') """

# adding the dividing line to the figure
#if (alpha==0.0):
plt.contour(X,Y,psi,levels=[0.0],colors='#CD2305',linewidths=2,linestyles='--')
plt.show()
