import matplotlib.pyplot as plt

from numpy import arange
from numpy.random import gamma
from math import cos
from math import sin
from math import atan
from math import sqrt
from random import random

#R0 SHAPE COMPONENTS FUNCTIONS
def a1(A20):
    return 4*pow(A20,3)/35

def a2(A20):
    return 6*pow(A20,2)/5

def K(A20):
    return a1(A20) + a2(A20) + 2

#R0 POSITION COMPONENTS FUNCTIONS
def b1(x,y,z):
    return pow(cos(atan(sqrt(x*x+y*y)/z)),2)
def b2(x,y,z,A20):
    return (1 + A20 *(0.5*(3*b1(x,y,z)-1)))

#GENERATE POSITIONS FUNCTION
def GenerateNucleonPosition(A20, R0, rhoBar, rhoRange, gamScale, gamShape):
    zConstraints = R0*(1+A20)
    sign = 1 if random() > 0.5 else -1
    z = gamma(gamShape,gamScale)
    while z >= zConstraints:
        z = gamma(gamShape, gamScale)
    x = min(rhoBar) + rhoRange*random()
    y = min(rhoBar) + rhoRange*random()
    R = R0 * b2(x,y,z,A20)
    R2 = R*R
    nucleonPosition = x*x + y*y + z*z
    loopCount = 1
    tries = 1
    while nucleonPosition > R2:
        x = min(rhoBar) + rhoRange*random()
        y = min(rhoBar) + rhoRange*random()
        R = R0 * b2(x,y,z,A20)
        R2 = R*R
        loopCount = loopCount + 1
        if loopCount > 10:
                z = gamma(gamShape,gamScale)
                while z >= zConstraints:
                    z = gamma(gamShape, gamScale)
                loopCount = 1
                tries = tries + 1
        if tries > 10000:
            return [x,y,sign*z,1]
    return [x,y,sign*z,0]

#CONSTANTS
PI = 3.14159265
Z = 92
N = 140
SCALE = 12

def Main():

    A = Z + N
    sampleSize = 1000
    iConst = 1
    A20 = 2
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
    shape = 1 + abs(max(zBar))/SCALE
    xProton = []
    zProton = []
    errorCount = 1
    for i in range(N):
        positions = GenerateNucleonPosition(A20,R0,rhoBar,rhoRange,SCALE,shape)
        xProton.append(positions[0])
        zProton.append(positions[2])
        if positions[3] != 0:
            errorCount = errorCount + 1
    print("Invalid positions: %i" %errorCount)
    plt.plot(zBar,rhoBar)
    plt.scatter(zProton,xProton,c='red')
    plt.show()


Main()