#PUBLIC LIBRARY IMPORTS
from math import sqrt
from numpy import arange
from numpy import poly1d
from numpy import polyfit
from statistics import mean

#USER DEFINED IMPORTS
from constants import PI
from constants import DATA_FOLDER_PATH
from utils import LoadJsonDataFile
from utils import SaveJsonDataFile
from utils import GetPositionsArray
from utils import FFunc
from utils import R0Func

#PRIVATE FUNCTIONS
def __GetCoulombEnergyFromPosition(Z,postProton):
    E = 0
    for i in range(Z):
        for j in range(i+1,Z):
            xDist = postProton[0][i] - postProton[0][j]
            yDist = postProton[1][i] - postProton[1][j]
            zDist = postProton[2][i] - postProton[2][j]
            r = sqrt(xDist*xDist + yDist*yDist + zDist*zDist)
            E = E + 1/r
    return E

def __GetSavedCoulombEnergies(N,Z,sampleCount):
    jsonKey = "coulomb_energy_array_{}_{}_{}".format(N+Z, Z, sampleCount)
    dataJson = LoadJsonDataFile()
    if jsonKey in dataJson:
        coulombEnergies = dataJson[jsonKey]
        return coulombEnergies
    return None

def __SaveCoulombEnergies(N, Z, sampleCount, coulombEnergies, overwrite=False):
    jsonKey = "coulomb_energy_array_{}_{}_{}".format(N+Z, Z, sampleCount)
    dataJson = LoadJsonDataFile()
    if overwrite or not jsonKey in dataJson:
        dataJson[jsonKey] = coulombEnergies
        SaveJsonDataFile(dataJson)

def __GetSavedFittedCoulombFunction(N,Z,sampleCount):
    jsonKey = "fitted_coeff_array_{}_{}_{}".format(N+Z, Z, sampleCount)
    dataJson = LoadJsonDataFile()
    if jsonKey in dataJson:
        savedCoeffList = dataJson[jsonKey]
        return poly1d(savedCoeffList)
    return None

def __SaveFittedCoulombFunction(N,Z,sampleCount, coeffList, overwrite = False):
    jsonKey = "fitted_coeff_array_{}_{}_{}".format(N+Z, Z, sampleCount)
    dataJson = LoadJsonDataFile()
    if overwrite or not jsonKey in dataJson:
        dataJson[jsonKey] = coeffList.tolist()
        SaveJsonDataFile(dataJson)

def __GetSavedSurfaceEnergies(A):
    jsonKey = "surface_energy_array_{}".format(A)
    dataJson = LoadJsonDataFile()
    if jsonKey in dataJson:
        return dataJson[jsonKey]
    return None

def __SaveSurfaceEnergies(A, surfaceEnergies, overwrite=False):
    jsonKey = "surface_energy_array_{}".format(A)
    dataJson = LoadJsonDataFile()
    if overwrite or not jsonKey in dataJson:
        dataJson[jsonKey] = surfaceEnergies
        SaveJsonDataFile(dataJson)

#ENERGY CALCULATION FUNCTIONS
def CalculateCoulombEnergies(N,Z,sampleCount):
    A = N + Z
    coulombEnergies = __GetSavedCoulombEnergies(N, Z, sampleCount)
    if coulombEnergies is not None:
        return coulombEnergies
    coulombEnergies, A20Plot = [], arange(0,2.01,0.01)
    for A20 in A20Plot:
        E = []
        for sample in range(sampleCount):
            postFile = open(DATA_FOLDER_PATH + r"{A}-{Z}/{sample}/{A}-{Z}_Proton_{A20:.2f}".format(A=A,Z=Z,A20=A20,sample=sample),'r')
            postProton = GetPositionsArray(postFile)
            E.append(__GetCoulombEnergyFromPosition(Z,postProton))
        coulombEnergies.append(mean(E))
    __SaveCoulombEnergies(N, Z, sampleCount, coulombEnergies)
    return coulombEnergies

def CalculateFittedCoulombEnergies(N,Z,sampleCount):
    coulombEnergies, A20Plot = CalculateCoulombEnergies(N,Z,sampleCount), arange(0,2.01,0.01)
    fitFunction = __GetSavedFittedCoulombFunction(N,Z,sampleCount)
    if fitFunction is None:
        fitFunction = poly1d(polyfit(A20Plot,coulombEnergies,9))
        __SaveFittedCoulombFunction(N, Z, sampleCount, fitFunction.c)
    fittedCoulombEnergies = []
    for A20 in A20Plot:
        fittedCoulombEnergies.append(fitFunction(A20))
    return fittedCoulombEnergies

def CalculateSurfaceEnergies(A):
    surfaceEnergies = __GetSavedSurfaceEnergies(A)
    if surfaceEnergies is None:
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
    __SaveSurfaceEnergies(A,surfaceEnergies)
    return surfaceEnergies

def CalculateTotalEnergies(N,Z,sampleCount):
    A = N+Z
    A_VOL = 15.75
    A_ASM = 23.7
    A_PAR = 11.18
    totalEnergies = []
    pairingTerm = 0

    if A % 2 != 1:
        if N % 2 == 0 and Z % 2 == 0:
            pairingTerm = A_PAR*pow(A,-0.5)
        else:
            pairingTerm = -A_PAR*pow(A,-0.5)
    
    volumeAsymmetryTerm = (A_VOL*A) - (A_ASM*pow((A-(2*Z)),2)/A)  

    surfaceEnergies = CalculateSurfaceEnergies(N+Z)
    coulombEnergies = CalculateFittedCoulombEnergies(N, Z, sampleCount)
    A_COL = (0.711*Z*Z)/(pow(N+Z,1/3)*coulombEnergies[0])
    A_SUR = 17.8/(4*PI*1.25*1.25)

    for i in range(len(surfaceEnergies)):
        totalEnergy = volumeAsymmetryTerm - (A_COL*coulombEnergies[i]) - (A_SUR*surfaceEnergies[i]) - pairingTerm
        totalEnergies.append(totalEnergy) 
    return totalEnergies

def CalculateFissionBarrier(N,Z,sampleCount):
    totalEnergies = CalculateTotalEnergies(N,Z,sampleCount)
    groundEnergy = totalEnergies[0]
    peakEnergy = max(totalEnergies)
    return peakEnergy - groundEnergy