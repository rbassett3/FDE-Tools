# FDE-Tools
Python code for computing fused density estimators (FDEs). Further details on fused density estimation can be found here. This package also contains a number auxillary functions for importing and building geometric networks from OpenStreetMaps (OSM) XML files.
## Installation & Configuration
1. Clone this repository to a local directory.
2. Install OSQP, and (optionally, for interior-point capabilities) CVXOPT
```
pip install OSQP CVXOPT
```
## Usage
### Univariate FDEs
Univariate FDEs can be declared with the syntax
```
fde = UnivarFDE((a,b), P)
```
where P is a collection of 1-D observations between a and b. For a penalty parameter lambda, the fde can be generated and solved via
```
fde.GenerateProblem()
fde.SolveProblem(lam)
```
Lastly, one can plot the fde with
```
fde.Plot()
plt.show()
```
Lastly, one can also perform cross-validation on the penalty parameter with 
```
fde.CrossValidate()
```
which returns a float.

### FDEs on Geometric Networks
Geometric network FDEs can be declared with the syntax
```
fde = FDE(L,P)
```
L is an array with elements as the segments of the geometric network. A segment is an array of 2D points. By default, these are assumed to be of the form (long, lat) and segment distances are calculated from these coordinates. As an example of this syntax,
```
L = [[[x1, y1], [x2, y2], [x3, y3]], [[x4, y4], [x2, y2], [x5, y5]]]
```
This example is a geometric network with two segments. The segments intersect at [x2, y2]. P is an array of points from L, so that 
```
P = [[x1, y1], [x3, y3], [x4, y4], [x5, y5]]
```
places observations at the endpoints of the above network. Generating, solving, plotting, and finding a lambda-parameter via cross-validation is all done as in the univariate case.
```
fde.GenerateProblem()
lam = fde.CrossValidate()
fde.SolveProblem(lam)
fde.Plot()
plt.show()
```

## Examples
### A univariate example
In this example, we will estimate a normal density with 100 data points.
```
import sys
sys.path.append("../../FDE-Tools") #import package
from FDE import *
import scipy.stats as stats

P = stats.norm.rvs(size = 100) #generate data
(a,b) = (-4,4)

fde = UnivarFDE((a,b), P) #Declare and solve FDE
fde.GenerateProblem()
fde.SolveProblem(.030)
fde.Plot()
plt.plot(x, stats.norm.pdf(x), color = 'darkorange')
plt.show()

#Perform cross-validation
print("Performing Cross Validation...")
lam = fde.CrossValidate()
fde.SolveProblem(lam)
print("Optimal lambda parameter: " + str(lam))
x = np.linspace(-4,4,100)
plt.plot(x, stats.norm.pdf(x), color = 'darkorange')
fde.Plot()
plt.show()
```
### A geometric network example
In this example we will load a geometric network from an OSM XML file ("SanDiego.xml"). We will also extract observations from the XML file as the locations of eateries in the downloaded region.
```
import sys
sys.path.append("../../FDE-Tools") #import package
from FDE import *
from HelperFunctions import *

#A window to focus on
lon = (-117.7, -117.149)
lat = (32.7071, 32.7216)

#Load data from xml, filter it to desired window
L = LoadMapXML("SanDiego.xml")
L = [FilterData(l, lon, lat) for l in L]
P = GetEateries("SanDiego.xml")
P = FilterData(P, lon, lat)

#Declare and solve the FDE
print("Declaring fused density estimator")
fde = FDE(L,P)
print("Generating Problem")
fde.GenerateProblem()
print("Solving Problem")
fde.SolveProblem(.022)
fde.Plot()
plt.show()
```          
