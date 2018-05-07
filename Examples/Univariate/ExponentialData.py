import sys
sys.path.append("../../FDE-Tools")
from FDE import *
import scipy.stats as stats

P = stats.expon.rvs(size = 100)
(a,b) = (0,6)
fde = UnivarFDE((a,b), P)
fde.GenerateProblem()
fde.SolveProblem(.03)

KeepGoing = input("Peform cross validation on this example (0/1)? ")
if KeepGoing == 1:
    print("Performing Cross Validation...")
    lam = fde.CrossValidate()
    fde.SolveProblem(lam)
    print("Optimal lambda parameter: " + str(lam))

fde.Plot()
x = np.linspace(-4,4,100)
plt.plot(x, stats.expon.pdf(x), color = 'darkorange')
plt.show()
