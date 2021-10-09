import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "..","..","FDE-Tools"))
from HelperFunctions import *
from FDE import *
import pandas as pd

#The latitude bounds in the Baghdad map are (33.2979, 33.3295).
#Longitude bounds are (44.3245, 44.3918).
lon = (44.3245, 44.3918)
lat = (33.2979, 33.3296)
print("Loading Map")
if os.path.isfile("BaghdadMap.npy"):
    L = np.load("BaghdadMap.npy",allow_pickle=True,encoding="latin1")
else:
    print (2)
    L = LoadMapXML('Baghdad.osm')
    L = [FilterData(l, lon, lat) for l in L] 
    np.save("BaghdadMap.npy", L)
print("Loading observations")
if os.path.isfile("GlobalTerrorismDatabase/BaghdadObs.npy"):
    P = np.load("GlobalTerrorismDatabase/BaghdadObs.npy")
else:
    D = pd.read_csv("GlobalTerrorismDatabase/gtd_13to16_0617dist.csv")
    IraqD = D[D.country == 95] #Country code for GTD data
    FilD = IraqD[np.logical_and(np.logical_and(IraqD.latitude >= lat[0], IraqD.latitude <= lat[1]), np.logical_and(IraqD.longitude <= lon[1], IraqD.longitude >= lon[0]))]
    P = np.array([FilD.longitude, FilD.latitude])
    P = P.T
    np.save("GlobalTerrorismDatabase/BigBaghdadObs.npy", P)
print("Declaring fused density estimator")
fde = FDE(L,P)
print("Generating Problem")
fde.GenerateProblem()
fde.Plot()
plt.savefig("Baghdad_JustObs.pdf")
print("Solving Problem")
fde.SolveProblem(.2)
fde.Plot()
plt.show()
