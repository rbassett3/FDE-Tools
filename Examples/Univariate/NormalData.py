import sys
sys.path.append("../../FDE-Tools")
from FDE import *
import scipy.stats as stats
random.seed(0)

P = stats.norm.rvs(size = 100)
(a,b) = (-4,4)
fde = UnivarFDE((a,b), P)
fde.GenerateProblem()
fde.SolveProblem(.030)
fde.Plot()

KeepGoing = eval(input("Perform cross validation on this example (1/0)? "))
if KeepGoing == 1:
    print("Performing Cross Validation...")
    lam = fde.CrossValidate()
    fde.SolveProblem(lam)
    print(("Optimal lambda parameter: " + str(lam)))
x = np.linspace(-4,4,100)
plt.plot(x, stats.norm.pdf(x), color = 'darkorange')
plt.show()
