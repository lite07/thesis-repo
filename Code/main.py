import multiprocessing as mp
import random as rnd
import time
import matplotlib.pyplot as plt
from generate_positions import GeneratePositions
from generate_graph import GenerateCoulombEnergyGraphs
from generate_graph import GenerateSurfaceEnergyGraphs
from calculate_energies import CalculateFissionBarrier

def main(N, Z, sampleCount):
    if len(N) < 2 and len(Z) <2:
        print("[Main] Please input a valid boundary")
    else:
        botN, botZ = N[0], Z[0]
        topN, topZ = N[1], Z[1]
        for neutronCount in range(botN, topN+1):
            for protonCount in range(botZ, topZ+1):
                print("[Main] Generating positions for {A}-{Z}.".format(A=neutronCount+protonCount,Z=protonCount))
                GeneratePositions(neutronCount, protonCount, sampleCount)
                #print("[Main] Generating coulomb graph for {A}-{Z} with {count} sample counts.".format(A=neutronCount+protonCount,Z=protonCount, count=sampleCount))
                #GenerateCoulombEnergyGraphs(neutronCount,protonCount,sampleCount)
                #print("[Main] Generating surface graph for {A}-{Z}.".format(A=neutronCount+protonCount,Z=protonCount))
                #GenerateSurfaceEnergyGraphs(neutronCount, protonCount)

#main([120,141],[80,100],15)

def test():
    NRange, ZRange = range(120,160), range(80,101)
    barrierPlot = []
    for N in NRange:
        barrierList = []
        for Z in ZRange:
            fissionBarrier = CalculateFissionBarrier(N,Z,10)
            barrierList.append(fissionBarrier)
            #print("{}-{}".format(N+Z,Z))
        barrierPlot.append(barrierList)
    plt.contourf(ZRange, NRange, barrierPlot)
    plt.colorbar()
    plt.xlabel("Proton (Z)")
    plt.ylabel("Neutron (N)")
    plt.title("Fission Barrier")
    plt.show()
    
#test()