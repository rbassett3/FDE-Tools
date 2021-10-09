import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "..","..","FDE-Tools"))
from HelperFunctions import *
from FDE import *
np.random.seed(0)

lon = (-121.9, -121.853)
lat = (36.587, 36.604)
print("Loading Map")
if os.path.isfile("Monterey.npy"):
    L = np.load("Monterey.npy",allow_pickle=True,encoding="latin1")
else:
    L = LoadMapXML("Monterey.xml")
    L = [FilterData(l, lon, lat) for l in L] 
    np.save("Monterey.npy", L)
print("Generating Observations")
P = (-121.88, 36.596) + (.007,.003)*np.random.randn(100,2)
print("Declaring fused density estimator")
fde = FDE(L,P)
print("Generating Problem")
fde.GenerateProblem()
print("Solving Problem")
fde.SolveProblem(.02)

#print("Performing Cross Validation...")
#lam = fde.CrossValidate()
#fde.SolveProblem(lam)
#print("Optimal lambda parameter: " + str(lam))

fde.Plot()
plt.show()
