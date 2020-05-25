import matplotlib.pyplot as plt
from numpy import array
from numpy import arange
from numpy import polyfit
from numpy import poly1d
from utils import CalculateCoulombEnergies

def GenerateCoulombEnergyGraphs(N,Z,sampleCount):
    EPlot, A20Plot = CalculateCoulombEnergies(N,Z,sampleCount), arange(0,2.01,0.01)
    fitFunction = poly1d(polyfit(A20Plot,EPlot,7))
    fittedCoulombEnergy = []
    for A20 in A20Plot:
        fittedCoulombEnergy.append(fitFunction(A20))
    plt.plot(A20Plot,fittedCoulombEnergy)
    plt.plot(A20Plot,EPlot)
    plt.savefig("Data/{A}-{Z}/Images/Coulomb_Graph.png".format(A=N+Z,Z=Z))

GenerateCoulombEnergyGraphs(143,92,15)
