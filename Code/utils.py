from constants import PI
from math import cos
from math import sin
from numpy import arange

#R0 SHAPE COMPONENTS FUNCTIONS
def a1(A20):
    return 4*pow(A20,3)/35

def a2(A20):
    return 6*pow(A20,2)/5

def K(A20):
    return a1(A20) + a2(A20) + 2

def GenerateCommonValues(A,A20):
    R0 = 1.25 * pow(2 * A / K(A20),1/3)
    dtheta = 0.00001*PI
    zBar, rhoBar, rBar = [], [], []
    for theta in arange(0,2*PI,dtheta):
        R = R0 * (1 + A20 * (0.5*(3*cos(theta)*cos(theta)-1)))
        zBar.append(R*cos(theta))
        rhoBar.append(R*sin(theta))
        rBar.append(R)

    return R0, zBar, rhoBar, rBar