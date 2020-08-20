from statistics import mean
import matplotlib.pyplot as plt
from numpy import arange
from utils import CalculateFittedCoulombEnergies
from utils import CalculateSurfaceEnergies
from generate_positions import GeneratePositions
from constants import PI
def CalculateTotalEnergy(N,Z):
    A = N + Z
    AV, AA = 14.6433, 21.0680
    
    energyVolAsymTerm = AV*A - AA*(A-2*Z)*(A-2*Z)/A   
    energyCoulomb = CalculateFittedCoulombEnergies(N,Z,45)
    energySurface = CalculateSurfaceEnergies(A)
    
    AS = 14.0788 / (4*PI*1.2257*1.2257)
    AC = 0.6442 * Z * Z / (pow(A,1/3)*energyCoulomb[0])
    totalEnergy = []
    for i in range(len(energyCoulomb)):
        totalEnergy.append(energyVolAsymTerm - (AS * energySurface[i] + AC * energyCoulomb[i]))
    A20Plot = arange(0,2.01,0.01)
    plt.plot(A20Plot,totalEnergy)
    plt.show()

CalculateTotalEnergy(143,92)

