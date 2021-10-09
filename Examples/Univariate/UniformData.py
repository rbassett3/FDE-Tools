import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "..","..","FDE-Tools"))
from FDE import *
import scipy.stats as stats

P = stats.uniform.rvs(size = 100)
(a,b) = (-0.1,1.1)
fde = UnivarFDE((a,b), P)
fde.GenerateProblem()
fde.SolveProblem(.03)

KeepGoing = input("Perform cross validation on this example (0/1)? ")
if KeepGoing == 1:
    print("Performing Cross Validation...")
    lam = fde.CrossValidate()
    FDE.SolveProblem(lam)
    print("Optimal lambda parameter: " + str(lam))

fde.Plot()
x = np.linspace(-0.1,1.1,100)
plt.plot(x, stats.uniform.pdf(x), color = 'darkorange')
plt.show()
