import matplotlib.pyplot as plt
from numpy import array
from numpy import arange
from numpy import polyfit
from numpy import poly1d
from utils import CalculateCoulombEnergies
from utils import CalculateSurfaceEnergies
from constants import DATA_FOLDER_PATH

def GenerateCoulombEnergyGraphs(N,Z,sampleCount):
    EPlot, A20Plot = CalculateCoulombEnergies(N,Z,sampleCount), arange(0,2.01,0.01)
    fitFunction = poly1d(polyfit(A20Plot,EPlot,7))
    fittedCoulombEnergy = []
    for A20 in A20Plot:
        fittedCoulombEnergy.append(fitFunction(A20))
    plt.plot(A20Plot,fittedCoulombEnergy)
    plt.plot(A20Plot,EPlot)
    plt.savefig(DATA_FOLDER_PATH + "{A}-{Z}/Images/Coulomb_Graph.png".format(A=N+Z,Z=Z))

#GenerateCoulombEnergyGraphs(143,92,45)

def GenerateSurfaceEnergyGraphs(N,Z):
    EPlot, A20Plot = CalculateSurfaceEnergies(N+Z), arange(0,2.01,0.01)
    fitFunction = poly1d(polyfit(A20Plot,EPlot,7))
    fittedSurfaceEnergy = []
    for A20 in A20Plot:
        fittedSurfaceEnergy.append(fitFunction(A20))
    plt.plot(A20Plot,fittedSurfaceEnergy)
    plt.plot(A20Plot,EPlot)
    plt.savefig(DATA_FOLDER_PATH + "{A}-{Z}/Images/Surface_Graph.png".format(A=N+Z,Z=Z))

GenerateSurfaceEnergyGraphs(143,92)