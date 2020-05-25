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
from os import makedirs
from random import random

#USER DEFINED IMPORTS
from constants import SCALE
from constants import PI
from utils import GenerateCommonValues
from utils import R0Func

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
def CreateFolder(A,Z,sampleCount=0,imageFolder = False):
    folderName = "Data/{A}-{Z}".format(A = A, Z = Z)
    if not path.exists(folderName):
        makedirs(folderName)
            
    if imageFolder :
        folderName = folderName + r"/Images"
        if not path.exists(folderName):
            makedirs(folderName)
    else:
        folderName = folderName + "/{sampleCount}".format(sampleCount = sampleCount)
        if not path.exists(folderName):
            makedirs(folderName)

def FileExist(filePath):
    if path.exists(filePath) : 
        return True
    return False

def OutputPositionsToFile(fileName, sampleNumber, postArray, nucleonCount):
    fileOutput = open(fileName,"w+")
    for i in range(nucleonCount):
        x, y, z = postArray[0][i],postArray[1][i],postArray[2][i]
        fileOutput.write("{x},{y},{z}\n".format(x=x,y=y,z=z))
    fileOutput.close()

def SavePositionFigure(fileName,zBar,rhoBar,postNeutron,postProton):
    plt.plot(zBar, rhoBar)
    plt.scatter(postNeutron[2],postNeutron[0],c='red',marker='.')
    plt.scatter(postProton[2],postProton[0],c='green',marker='.')
    plt.savefig(fileName)
    plt.clf()

def GeneratePositions(N,Z,sampleCount):
    A = Z + N 
    for sample in range(sampleCount):
        CreateFolder(A,Z,sample)
        CreateFolder(A,Z,imageFolder = True)
        for A20 in arange(0,2.01,0.01):
            neutronFileName = "Data/{A}-{Z}/{sampleNumber}/{A}-{Z}_{type}_{A20:.2f}".format(A = A, Z = Z, sampleNumber = sample, type = "Neutron", A20=A20)
            protonFileName = "Data/{A}-{Z}/{sampleNumber}/{A}-{Z}_{type}_{A20:.2f}".format(A = A, Z = Z, sampleNumber = sample, type = "Proton", A20=A20)
            imageFileName = "Data/{A}-{Z}/Images/{A20:.2f}.png".format(A = A, Z = Z, A20 = A20)
            if FileExist(neutronFileName) and FileExist(imageFileName) and FileExist(protonFileName):
                continue
            
            R0 = R0Func(A20,A)
            zBar, rhoBar, rBar = GenerateCommonValues(A,A20)
            rhoRange = max(rhoBar) - min(rhoBar)
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

            OutputPositionsToFile(neutronFileName,sample,postNeutron,A-Z)
            OutputPositionsToFile(protonFileName,sample,postProton,Z)
            if sample == 0:
                SavePositionFigure(imageFileName,zBar,rhoBar,postNeutron,postProton)
    

N = 143 #int(sys.argv[1])
Z = 92 #int(sys.argv[2])
sampleCount = 15
GeneratePositions(N,Z,sampleCount)