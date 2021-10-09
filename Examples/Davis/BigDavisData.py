import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "..","..","FDE-Tools"))
from FDE import *
from HelperFunctions import *

lat = (-121.795, -121.674)
lon = (38.5282, 38.575)
print("This is a large scale example. Initially loading the map data can take a while.")
print("Loading Map")
if os.path.isfile("BigDavisMap.npy"):
    L = np.load("BigDavisMap.npy",allow_pickle=True,encoding="latin1")
else:
    L = LoadMapXML("BigDavisMap.xml")
    L = [FilterData(l, lat, lon) for l in L]
    np.save("BigDavisMap.npy", L)
    
print("Loading locations of eateries")
if os.path.isfile("BigDavisEateries.npy"):
    P = np.load("BigDavisEateries.npy")
else:
    P = GetEateries("BigDavisMap.xml")
    P = FilterData(P, lat, lon) 
    np.save("BigDavisEateries.npy", P)

print("Defining fused density estimator")
fde = FDE(L,P)
print("Generating Problem")
fde.GenerateProblem()
print("Solving Problem")
fde.SolveProblem(.0271)
fde.Plot()
plt.show()
