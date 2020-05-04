import matplotlib.pyplot as plt
from numpy import arange
from math import cos
from math import sin

PI = 3.14159265
Z = 92
N = 140
def a1(A20):
    return 4*pow(A20,3)/35

def a2(A20):
    return 6*pow(A20,2)/5

def K(A20):
    return a1(A20) + a2(A20) + 2

def Main():

    A = Z + N
    sampleSize = 1000
    iConst = 1
    A20 = 0
    R0 = 1.2 * pow(2 * A / K(A20),1/3)
    dx = 0.001
    dy = 0.001
    dz = 0.001
    dtheta = 0.001*PI
    dphi = 0.001*PI
    zBar = []
    rhoBar = []
    rBar = []
    for theta in arange(0,2*PI,dtheta):
        R = R0 * (1 + A20 * (0.5*(3*cos(theta)*cos(theta)-1)))
        zBar.append(R*cos(theta))
        rhoBar.append(R*sin(theta))
        rBar.append(R)
    zRange = max(zBar) - min(zBar)
    rhoRange = max(rhoBar) - min(rhoBar)
    rRange = max(rBar) - min(rBar)
    plt.plot(zBar,rhoBar)
    plt.show()

Main()