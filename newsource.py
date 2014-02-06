import numpy as np
import matplotlib.pyplot as plt
from math import *
N = 200                            # Number of points in each direction
xStart,xEnd = -5.0,5.0            # x-direction boundaries
yStart,yEnd = -5.0,5.0            # y-direction boundaries
x = np.linspace(xStart,xEnd,N)    # x 1D-array
y = np.linspace(yStart,yEnd,N)    # y 1D-array

[X,Y] = np.meshgrid(x,y)            # generation of the mesh grid

dx=(xEnd-xStart)/N
dy=(yEnd-yStart)/N

# plotting the mesh
size = 10
"""
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.scatter(X,Y,s=10,c='#CD2305',marker='o',linewidth=0.1)
plt.show()
"""
strengthSource = 5.0                   # source strength
strengthSink=-5.0
xSource, ySource = -0.5,0.0            # location of the source
xSink, ySink=.5,0.0

uSource = np.empty((N,N),dtype=float)  # creation of a 2D-array for u
vSource = np.empty((N,N),dtype=float)  # creation of a 2D-array for v
uSink = np.empty((N,N),dtype=float)  # creation of a 2D-array for u
vSink = np.empty((N,N),dtype=float)  # creation of a 2D-array for v..
u = np.empty((N,N),dtype=float) 
v= np.empty((N,N),dtype=float) 
Res= np.empty((N,N),dtype=float) 

U_inf=1.0
V_inf=0.0

phisource=np.empty((N,N),dtype=float)
phisink=np.empty((N,N),dtype=float)
phifreestream=np.empty((N,N),dtype=float)
phi=np.empty((N,N),dtype=float)



# computing the velocity components at every point on the mesh grid

uSource=(strengthSource/(2*pi))*((X-xSource)/((X-xSource)**2+(Y-ySource)**2))

vSource= (strengthSource/(2*pi))*((Y-ySource)/((X-xSource)**2+(Y-ySource)**2))

uSink=(strengthSink/(2*pi))*((X-xSink)/((X-xSink)**2+(Y-ySink)**2))

vSink= (strengthSink/(2*pi))*((Y-ySink)/((X-xSink)**2+(Y-ySink)**2))

u=uSource+uSink+U_inf
v=vSource+vSink+V_inf

for i in range(N):
    for j in range(N):
        phisource[i][j]=(strengthSource/(2*pi))*log((X[i][j]-xSource)**2+(Y[i][j]-ySource)**2)
        phisink[i][j]=(strengthSink/(2*pi))*log((X[i][j]-xSink)**2+(Y[i][j]-ySink)**2)
        phifreestream[i][j]=U_inf*X[i][j]

phi=phisource+phisink+phifreestream

# plotting the streamlines
size = 10

for i in range(N):
    for j in range(N):
        Res[i][j]=sqrt(u[i][j]**2+v[i][j]**2)

plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=2,arrowstyle='->')
plt.scatter(xSource,ySource,c='#CD2305',s=80,marker='o')
plt.scatter(xSink,ySink,c='#CD2305',s=80,marker='o')
plt.show()

Cp = 1.0-(u**2+v**2)/U_inf**2

plt.figure()
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
contf = plt.contour(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
"""cbar = plt.colorbar(contf)
cbar.set_label(r'$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0]
plt.scatter(xSource,ySource,c='#CD2305',s=80,marker='o')
plt.scatter(xSink,ySink,c='#CD2305',s=80,marker='o')"""
plt.scatter([xSource,xSink],[ySource,ySink],c='#CD2305',s=80,marker='o')
plt.contour(X,Y,phi,levels=[0.0],colors='#CD2305',linewidths=2,linestyles='solid')
plt.show()


