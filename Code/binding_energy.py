import constants
import csv
import matplotlib.pyplot as pyPlot

def CalculateBE(A, Z):
    Ev = constants.VOL_CONSTANT*A
    Es = constants.SFC_CONSTANT*(A**(2/3))
    Ec = constants.CLB_CONSTANT*Z*(Z-1)*(A**(-1/3))
    Ea = constants.ASM_CONSTANT*((A-2*Z)**2)/A
    Ep = 0
    if A%2 == 1:
        Ep = 0
    elif Z%2 == 1:
        Ep = -constants.PAI_CONSTANT*(A**(-3/4))   
    else:
        Ep = constants.PAI_CONSTANT*(A**(-3/4))  
    return Ev + Es + Ec + Ea + Ep

def GetElements():
    result = list()
    with open('elements_list.csv') as elementsFile:
        elements = csv.reader(elementsFile,delimiter=',')
        topRow = True
        for element in elements:
            if not topRow:
                toAdd = [int(element[4]) + int(element[5]),int(element[5])]
                result.append(toAdd)
            else:
                topRow = False
    return result

energyPerMass = list()
massNumber = list()

for element in GetElements():
    massNumber.append(element[0])
    energyPerMass.append(-CalculateBE(element[0],element[1])/element[0])

pyPlot.plot(massNumber, energyPerMass)
pyPlot.xlabel("Mass Number")
pyPlot.ylabel("Binding Energy per Mass (MeV)")
pyPlot.grid(all)
pyPlot.show()