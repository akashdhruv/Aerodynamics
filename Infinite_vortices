import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 100 
xStart,xEnd = -4.0,4.0
yStart,yEnd = -2.0,2.0 
x = np.linspace(xStart,xEnd,N)
y = np.linspace(yStart,yEnd,N)    
X,Y = np.meshgrid(x,y)          

Uinf = 1.0
alphaInDegrees = 0.0
alpha = alphaInDegrees*pi/180


uFreestream = Uinf*cos(alpha)
vFreestream = Uinf*sin(alpha)


psiFreestream = + Uinf*cos(alpha)*Y - Uinf*sin(alpha)*X

def flowvortexinfinite(ga,a,xn,yn):
    u = (+ga/(2*pi))*((sinh(2*pi*yn/a))/((cosh(2*pi*yn/a)-cos(2*pi*xn/a))))
    v = (-ga/(2*pi))*((sin(2*pi*yn/a))/((cosh(2*pi*yn/a)-cos(2*pi*xn/a))))
    return u,v
    
a=0.41
ga=.001

uvortex=np.zeros((N,N),dtype=float)
vvortex=np.zeros((N,N),dtype=float)
psivortex=np.zeros((N,N),dtype=float)

for i in range(N):
    for j in range(N):
        uvortex[i][j],vvortex[i][j]=flowvortexinfinite(ga,a,X[i][j],Y[i][j])

    


size = 10
plt.figure()
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uvortex,vvortex,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.show()  