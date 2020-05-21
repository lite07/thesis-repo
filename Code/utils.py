from constants import PI
from math import cos
from math import sin
from math import sqrt
from numpy import arange

#PRIVATE FUNCTIONS
def __a1(A20):
    return 4*pow(A20,3)/35

def __a2(A20):
    return 6*pow(A20,2)/5

def __K(A20):
    return __a1(A20) + __a2(A20) + 2

def __ConvertToDouble(postString):
    splitter = postString.split(',')
    return float(splitter[0]), float(splitter[1]) , float(splitter[2])

#PUBLIC UTILITIES FUNCTION
def GenerateCommonValues(A,A20):
    R0 = 1.25 * pow(2 * A / __K(A20),1/3)
    dtheta = 0.00001*PI
    zBar, rhoBar, rBar = [], [], []
    for theta in arange(0,2*PI,dtheta):
        R = R0 * (1 + A20 * (0.5*(3*cos(theta)*cos(theta)-1)))
        zBar.append(R*cos(theta))
        rhoBar.append(R*sin(theta))
        rBar.append(R)

    return R0, zBar, rhoBar, rBar

def GetPositionsArray(postFile):
    Positions = [[],[],[]]
    for postString in postFile:
        x, y, z = __ConvertToDouble(postString)
        Positions[0].append(x)
        Positions[1].append(y)
        Positions[2].append(z)
    return Positions

def GetCoulombEnergyFromPosition(Z,postProton):
    E = 0
    for i in range(Z):
        for j in range(i+1,Z):
            xDist = postProton[0][i] - postProton[0][j]
            yDist = postProton[1][i] - postProton[1][j]
            zDist = postProton[2][i] - postProton[2][j]
            r = sqrt(xDist*xDist + yDist*yDist + zDist*zDist)
            E = E + 1/r
    return E