#PUBLIC LIBRARY IMPORTS
import matplotlib.pyplot as plt
from numpy import array
from numpy import arange
from numpy import polyfit
from numpy import poly1d

#USER DEFINED IMPORTS
from calculate_energies import CalculateCoulombEnergies
from calculate_energies import CalculateSurfaceEnergies
from utils import FileExist
from constants import DATA_FOLDER_PATH

def GenerateCoulombEnergyGraphs(N,Z,sampleCount):
    fileName = DATA_FOLDER_PATH + "{A}-{Z}/Images/Coulomb_Graph_{count}.png".format(A=N+Z,Z=Z,count=sampleCount)
    if not FileExist(fileName):
        EPlot, A20Plot = CalculateCoulombEnergies(N,Z,sampleCount), arange(0,2.01,0.01)
        fitFunction = poly1d(polyfit(A20Plot,EPlot,7))
        fittedCoulombEnergy = []
        for A20 in A20Plot:
            fittedCoulombEnergy.append(fitFunction(A20))
        plt.plot(A20Plot,fittedCoulombEnergy)
        plt.plot(A20Plot,EPlot)
        plt.savefig(fileName)
        plt.clf()

def GenerateSurfaceEnergyGraphs(N,Z):
    fileName = DATA_FOLDER_PATH + "{A}-{Z}/Images/Surface_Graph.png".format(A=N+Z,Z=Z)
    if not FileExist(fileName):
        EPlot, A20Plot = CalculateSurfaceEnergies(N+Z), arange(0,2.01,0.01)
        fitFunction = poly1d(polyfit(A20Plot,EPlot,7))
        fittedSurfaceEnergy = []
        for A20 in A20Plot:
            fittedSurfaceEnergy.append(fitFunction(A20))
        plt.plot(A20Plot,fittedSurfaceEnergy)
        plt.plot(A20Plot,EPlot)
        plt.savefig(fileName)
        plt.clf()

GenerateCoulombEnergyGraphs(143,92,10)