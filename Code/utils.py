#PUBLIC LIBRARY IMPORTS
import matplotlib.pyplot as plt
import json
from array import *
from math import cos
from math import sin
from math import sqrt
from numpy import arange
from numpy import polyfit
from numpy import poly1d
from os import path
from statistics import mean

#USER DEFINED IMPORTS
from constants import PI
from constants import DATA_FOLDER_PATH

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
def FileExist(filePath):
    if path.exists(filePath) : 
        return True
    return False

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

def LoadJsonDataFile():
    jsonDataFileName = DATA_FOLDER_PATH + "data.json"
    if not FileExist(jsonDataFileName):
        print("[Utils.LoadJsonDataFile] 'data.json' is not found. Creating.")
        emptyJson = {}
        with open(jsonDataFileName, 'w') as jsonDataFile:
            json.dump(emptyJson, jsonDataFile)
        return emptyJson
    with open(jsonDataFileName) as jsonDataFile:
        dataJson = json.load(jsonDataFile)
        return dataJson

def SaveJsonDataFile(dataJson):
    jsonDataFileName = DATA_FOLDER_PATH + "data.json"
    with open(jsonDataFileName, 'w') as jsonDataFile:
        json.dump(dataJson, jsonDataFile, indent=4, sort_keys=True)