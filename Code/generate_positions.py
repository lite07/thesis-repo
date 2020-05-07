#PACKAGE IMPORTS
import matplotlib.pyplot as plt
import sys
from numpy import arange
from numpy.random import gamma
from math import cos
from math import sin
from math import atan
from math import sqrt
from os import path
from os import mkdir
from random import random

#USER DEFINED IMPORTS
from constants import SCALE
from constants import PI
from utils import GenerateCommonValues

#R0 POSITION COMPONENTS FUNCTIONS
def b1(x,y,z):
    return pow(cos(atan(sqrt(x*x+y*y)/z)),2)
def b2(x,y,z,A20):
    return (1 + A20 *(0.5*(3*b1(x,y,z)-1)))

#GENERATE POSITIONS FUNCTION
def GenerateNucleonPosition(A20, R0, rhoMin, rhoRange, gamScale, gamShape):
    zConstraints = R0*(1+A20)
    sign = 1 if random() > 0.5 else -1
    z = gamma(gamShape,gamScale)
    while z >= zConstraints:
        z = gamma(gamShape, gamScale)
    x = rhoMin + rhoRange*random()
    y = rhoMin + rhoRange*random()
    R = R0 * b2(x,y,z,A20)
    loopCount = 1
    while x*x + y*y + z*z >= R*R:
        x = rhoMin + rhoRange*random()
        y = rhoMin + rhoRange*random()
        R = R0 * b2(x,y,z,A20)
        loopCount = loopCount + 1
        if loopCount > 5000:
                z = gamma(gamShape,gamScale)
                while z >= zConstraints:
                    z = gamma(gamShape, gamScale)
                loopCount = 1
    return [x,y,sign*z]

#FOLDERING AND OUTPUT
def CreateFolder(A,Z,imageFolder = False):
    folderName = "{A}-{Z}".format(A = A, Z = Z)
    if not path.exists(folderName):
        mkdir(folderName)
    if imageFolder :
        folderName = folderName + r"/Images"
        if not path.exists(folderName):
            mkdir(folderName)

def OutputPositionsToFile(A20, A, Z, postArray, nucleonType):
    CreateFolder(A,Z)
    fileOutput = open("{A}-{Z}/{A}-{Z}_{type}_{A20:.2f}".format(A = A, Z = Z, type = nucleonType, A20=A20),"w+")
    nucleonCount = Z if nucleonType == "Proton" else A-Z
    for i in range(nucleonCount):
        x, y, z = postArray[0][i],postArray[1][i],postArray[2][i]
        fileOutput.write("{x},{y},{z}\n".format(x=x,y=y,z=z))
    fileOutput.close()

def SavePositionFigure(A20,A,Z,zBar,rhoBar,postNeutron,postProton):
    CreateFolder(A,Z,imageFolder=True)
    plt.plot(zBar, rhoBar)
    plt.scatter(postNeutron[2],postNeutron[0],c='red',marker='.')
    plt.scatter(postProton[2],postProton[0],c='green',marker='.')
    plt.savefig("{A}-{Z}/Images/{A20:.2f}.png".format(A = A, Z = Z, A20 = A20))
    plt.clf()

#CONSTANTS
Z = 92
N = 140

def GeneratePositions(N,Z):
    A = Z + N 
    for A20 in arange(0,2.01,0.01):
        R0, zBar, rhoBar, rBar = GenerateCommonValues(A,A20)
        rhoRange = max(rhoBar) - min(rhoBar)
        rRange = max(rBar) - min(rBar)
        shape = 1 + abs(max(zBar))/SCALE

        postNeutron = [[],[],[]]
        for i in range(N):
            newPositions = GenerateNucleonPosition(A20,R0,min(rhoBar),rhoRange,SCALE,shape)
            postNeutron[0].append(newPositions[0])
            postNeutron[1].append(newPositions[1])
            postNeutron[2].append(newPositions[2])

        postProton = [[],[],[]]
        for i in range(Z):
            newPositions = GenerateNucleonPosition(A20,R0,min(rhoBar),rhoRange,SCALE,shape)
            postProton[0].append(newPositions[0])
            postProton[1].append(newPositions[1])
            postProton[2].append(newPositions[2])

        OutputPositionsToFile(A20,A,Z,postNeutron,"Neutron")
        OutputPositionsToFile(A20,A,Z,postProton,"Proton")
        SavePositionFigure(A20,A,Z,zBar,rhoBar,postNeutron,postProton)

N = int(sys.argv[1])
Z = int(sys.argv[2])

GeneratePositions(N,Z)