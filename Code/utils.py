from constants import PI
from constants import DATA_FOLDER_PATH
from math import cos
from math import sin
from math import sqrt
from numpy import arange
from numpy import polyfit
from numpy import poly1d
from statistics import mean

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
def R0Func(A20, A):
    return 1.25 * pow(2 * A / __K(A20),1/3)

def PFunc(theta):
    return 0.5*(3*cos(theta)*cos(theta) - 1)

def RFunc(R0, A20, theta):
    return  R0 * (1 + A20 * PFunc(theta))

def DerivPFunc(theta):
    return 3*cos(theta)

def FFunc(R0, A20, theta):
    R = RFunc(R0,A20,theta)
    return R*sin(theta)*sqrt(pow(R0*A20*DerivPFunc(theta)*sin(theta),2) + pow(R,2)) 

def GenerateCommonValues(A,A20):
    R0 = R0Func(A20, A)
    dtheta = 0.0001*PI
    zBar, rhoBar, rBar = [], [], []
    for theta in arange(0,2*PI + dtheta,dtheta):
        R = RFunc(R0, A20, theta)
        zBar.append(R*cos(theta))
        rhoBar.append(R*sin(theta))
        rBar.append(R)

    return zBar, rhoBar, rBar

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

def CalculateCoulombEnergies(N,Z,sampleCount):
    A = N + Z
    coulombEnergies, A20Plot = [], arange(0,2.01,0.01)
    for A20 in A20Plot:
        E = []
        for sample in range(sampleCount):
            postFile = open(DATA_FOLDER_PATH + r"{A}-{Z}/{sample}/{A}-{Z}_Proton_{A20:.2f}".format(A=A,Z=Z,A20=A20,sample=sample),'r')
            postProton = GetPositionsArray(postFile)
            E.append(GetCoulombEnergyFromPosition(Z,postProton))
        coulombEnergies.append(mean(E))
    return coulombEnergies

def CalculateFittedCoulombEnergies(N,Z,sampleCount):
    A = N + Z
    coulombEnergies, A20Plot = CalculateCoulombEnergies(N,Z,sampleCount), arange(0,2.01,0.01)
    fitFunction = poly1d(polyfit(A20Plot,coulombEnergies,7))
    fittedCoulombEnergies = []
    for A20 in A20Plot:
        fittedCoulombEnergies.append(fitFunction(A20))
    return fittedCoulombEnergies

def CalculateSurfaceEnergies(A):
    nPartition = 1000
    dTheta = PI/nPartition
    A20Plot = arange(0,2.01,0.01)
    surfaceEnergies = []
    for A20 in A20Plot:
        R0 = R0Func(A20,A)
        surfaceEnergy = 0
        theta = 0
        for i in range(nPartition):
            if i == 0 :
                surfaceEnergy = surfaceEnergy + FFunc(R0,A20,theta)
            elif i == nPartition-1:
                surfaceEnergy = surfaceEnergy + FFunc(R0,A20,theta)
            elif i % 2 == 0:
                surfaceEnergy = surfaceEnergy + 2 * FFunc(R0,A20,theta)
            else:
                surfaceEnergy = surfaceEnergy + 4 * FFunc(R0,A20,theta)
            theta = theta + dTheta
        surfaceEnergies.append(2*PI*dTheta*surfaceEnergy*(1/3))
    return surfaceEnergies