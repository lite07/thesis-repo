import matplotlib.pyplot as plt
from numpy import array
from numpy import arange
from numpy import polyfit
from numpy import poly1d
from math import sqrt
from statistics import mean
from utils import GetPositionsArray
from utils import GetCoulombEnergyFromPosition

def GenerateCoulombEnergyGraphs(N,Z,sampleCount):
    A = N + Z
    EPlot, A20Plot = [], []
    for A20 in arange(0,2.01,0.01):
        E = []
        for sample in range(sampleCount):
            if sample == 0 : 
                A20Plot.append(A20)
            postFile = open(r"Data/{A}-{Z}/{sample}/{A}-{Z}_Proton_{A20:.2f}".format(A=A,Z=Z,A20=A20,sample=sample),'r')
            postProton = GetPositionsArray(postFile)
            E.append(GetCoulombEnergyFromPosition(Z,postProton))
        EPlot.append(mean(E))
    fitFunction = poly1d(polyfit(A20Plot,EPlot,7))
    fittedCoulombEnergy = []
    for A20 in A20Plot:
        fittedCoulombEnergy.append(fitFunction(A20))
    plt.plot(A20Plot,fittedCoulombEnergy)
    plt.plot(A20Plot,EPlot)
    plt.savefig("Data/{A}-{Z}/Images/Coulomb_Graph.png".format(A=A,Z=Z))

GenerateCoulombEnergyGraphs(142,90,15)
