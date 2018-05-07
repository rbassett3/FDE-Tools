import sys
sys.path.append("../../FDE-Tools")
from FDE import *
from HelperFunctions import *

lat = (-121.74845, -121.73161)
lon = (38.54081, 38.54824)

print("Loading Map")
if os.path.isfile("DavisMap.npy"):
    L = np.load("DavisMap.npy")
else:
    L = LoadMapXML("DavisMap.osm")
    L = [FilterData(l, lat, lon) for l in L] 
    np.save("DavisMap.npy", L)

print("Loading locations of eateries")
if os.path.isfile("DavisEateries.npy"):
    P = np.load("DavisEateries.npy")
else:
    P = GetEateries("DavisMap.osm")
    P = FilterData(P, lat, lon)
    np.save("DavisEateries.npy", P)

print("Defining fused density estimator")
fde = FDE(L,P)
print("Generating Problem")
fde.GenerateProblem()
print("Solving Problem")
fde.SolveProblem(.11)
fde.Plot()



