import pdb
import numpy as np
import scipy.stats
import scipy.sparse as sp
import scipy.sparse.csgraph as csgraph
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import collections
import math
from timeit import default_timer as timer
import random
import osqp
try:
    import osqp
except ImportError:
    ImportError("Default solver OSQP is not installed. Please install before continuing")
try:
    import cvxopt #An alternative solver for dual
except ImportError:
    pass
try:
    import cvxpy as cvx #Required to solve primal problem.
except ImportError: 
    pass

class FDE: 
    def __init__(self, Lines, Obs, OneConComp = True, lat_lon = True):
        self.lat_lon = lat_lon
        #Map to tuples for use with dictionary
        Ltup = [list(map(tuple, l)) for l in Lines]
        self.Ltup = Ltup
        #A collection of all points in L 
        LumpL = np.vstack(Lines) #set -> list removes duplicates
        LumpL = list(map(tuple, LumpL)) 
        LumpL = set(LumpL)
        LumpL = list(LumpL)
        #Dictionary connecting points to the corresponding node
        self.Point2Node = dict(list(zip(list(map(tuple, LumpL)), list(range(len(LumpL))))))
        #Build Adjacency matrix, and dictionaries connecting segments to length and nodes to length
        (Adj, self.Seg2Node, self.Node2Len) = self.BuildAdjMat(Ltup, self.Point2Node)
        #In Compressed sparse column format because COO format does not permit column slicing
        self.AdjCsc = Adj.tocsc()
        if OneConComp == True:
            (Adj, self.Point2Node, self.Seg2Node, self.Node2Len) = self.BigConnComp(self.AdjCsc, self.Point2Node, self.Seg2Node, self.Node2Len)
            #Form new LumpL and project onto it
            LumpL = list(self.Point2Node.keys())
            list(set(LumpL))
            self.AdjCsc = Adj.tocsc()
        Obs = self.ProjectToNetwork(Obs, LumpL) #The projection
    
        #Updating the dictionaries
        self.Point2Obs = dict(list(zip(LumpL, np.zeros(len(LumpL)))))
        self.Point2Obs.update(collections.Counter([tuple(obs) for obs in Obs]))
        self.numNodes = self.AdjCsc.shape[0]
        Node2Point = dict(list(zip(list(self.Point2Node.values()), list(self.Point2Node.keys()))))
        #If node is not a point, it is [False] in Node2Point dictionary
        Node2Point.update(dict(list(zip(list(range(len(self.Point2Node), self.numNodes)), [False]*(self.numNodes-len(LumpL))))))
        #False has no observations, so Node2Point[Point2Obs[n]] = 0 when n corresponds to a segment
        self.Point2Obs[False] = 0
        #A dictionary from nodes to observations
        self.Node2Obs = dict(list(zip(list(range(self.numNodes)), [self.Point2Obs[Node2Point[n]] for n in range(self.numNodes)])))
    
        self.Problem_Declared = False
        self.Problem_Solved = False

    def GenerateProblem(self):
        #Find the points which are important--those which are an observation or have degree >= 3
        ImpNodes = [(self.Node2Obs[n] != 0) or (self.AdjCsc[:, n].nnz >= 3) for n in range(self.numNodes)]
        numImp = sum(ImpNodes)
    
        #Remove those points from the graph
        RemI = sp.diags(np.logical_not(ImpNodes), dtype = np.bool)
        #Find the connect components of the graph with with the extraneous nodes removed.
        #We remove these nodes because theory tell us that no break point can occur there.
        (numComp, ConComp) = csgraph.connected_components(RemI.dot(self.AdjCsc.dot(RemI)))
    
        #Construct a quotient graph, where we group the unimportant nodes into a single node.
        #Linear in number of edges. Thanks Kirill!
        Edges = list(zip(self.AdjCsc.nonzero()[0], self.AdjCsc.nonzero()[1]))
        NewEdges = np.array([[ConComp[e[0]], ConComp[e[1]]] for e in Edges if ConComp[e[0]] != ConComp[e[1]]])
        #Reduced adjacency matrix
        RedAdj = sp.coo_matrix(([True]*NewEdges.shape[0], (NewEdges[:,0], NewEdges[:,1])))
    
        #Build the oriented edge incidence matrix
        D = self.BuildEdgeIncidence(RedAdj).tocsc()
    
        #A matrix of booleans. Rows are connected components, columns are nodes. 
        #True if node is in connected component.
        CompMat = sp.csr_matrix([ConComp == ind for ind in range(numComp)])
        #Compute the length of each connected component by summing over nodes in component.
        LenVec = np.array([self.Node2Len[i] for i in range(len(self.Node2Len))])
        #The number of observations of each connected component
        ObsVec = np.array([self.Node2Obs[i] for i in range(len(self.Node2Obs))])
        #Updating length and observation dictionaries to the quotient.
        RedNode2Len = dict(list(zip(list(range(numComp)), [(LenVec[CompMat[ind, :].nonzero()[1]]).sum() for ind in range(numComp)])))
        RedNode2Obs = dict(list(zip(list(range(numComp)), [(ObsVec[CompMat[ind, :].nonzero()[1]]).sum() for ind in range(numComp)])))
    
        #Vectors of these quantities
        RedLenVec = np.array(list(RedNode2Len.values()))
        RedObsVec = np.array(list(RedNode2Obs.values()))
    
        s_inds = RedLenVec != 0
        u_inds = np.logical_not(s_inds)
        self.D1 = D[:, s_inds]
        self.D2 = D[:, u_inds]
        self.s = RedLenVec[s_inds]
        self.u = -RedObsVec[u_inds]/(RedObsVec[u_inds].sum())
    
        #Seg -> OldNode -> NewNode -> Index among columns of D (same as indices of s)
        OldNode2NewNode = dict(list(zip(list(range(self.numNodes)), [ConComp[i] for i in range(self.numNodes)])))
        NewNode2s_ind = dict(list(zip(s_inds.nonzero()[0], list(range(len(self.s))))))
        self.Seg2s_index = dict(list(zip(list(self.Seg2Node.keys()), [NewNode2s_ind[OldNode2NewNode[self.Seg2Node[seg]]] for seg in list(self.Seg2Node.keys())]))) 
        self.Problem_Declared = True
    
    def ProjectToNetwork(self, P, LumpL):
        '''ProjectToNetwork projects each point in P onto its closest point in LumpL''' 
        NewP = [LumpL[np.argmin(np.linalg.norm(LumpL - p, ord = 2, axis = 1))] for p in P]
        #Could use map instead of list comprehension
        #NewP = map(lambda p: LumpL[np.argmin(np.linalg.norm(LumpL - p, ord = 2, axis = 1))], P)
        return NewP
    
    def BuildAdjMat(self, Ltup, Point2Node):
        '''Builds the geometric network (via adjacency matrix). There are many irrelevant points in here. 
        
        Both points and segments are represented as nodes in this adjacency matrix.
        Hence a graph of of the form:     o-----o      has 3 nodes. 
        1 -> left point. 2 -> edge. 3 -> right point.
        The adjacency matrix is [[0, 1, 0], [1, 0, 1], [0, 1, 0]].
        
        BuildAdjMat returns a dictionary of segments to nodes, an array of nodes to lengths (is 0 when the node is a point), and a 0-1 adjacency matrix.'''
        i = []
        j = []
        data = []
        lengths = []
        numSegs = sum([len(l)-1 for l in Ltup])
        Segs = []
        list(map(Segs.extend, [list(zip(l[:-1], l[1:])) for l in Ltup]))
        Seg2Node = dict(list(zip(Segs, list(range(len(Point2Node), len(Point2Node)+len(Segs))))))
        if self.lat_lon == False: #treat the tuples as elements in R^n
            #pdb.set_trace()
            Node2Len = dict(list(zip(list(range(len(Point2Node)+len(Segs))), [0]*len(Point2Node)+[np.linalg.norm(np.subtract(s[1],s[0])) for s in Segs])))
        else: 
            #treat the tuples as (lon,lat) coordinates 
            #Use Haversine formula to compute distance
            def haversine(x1, x2, miles = True):
                '''The Haversine formula. This function is based on the 
                haversine python package.'''
                AVG_EARTH_RADIUS = 6371  # in km
                MILES_PER_KILOMETER = 0.621371

                # unpack latitude/longitude
                lng1, lat1 = x1
                lng2, lat2 = x2
        
                # convert all latitudes/longitudes from decimal degrees to radians
                lat1, lng1, lat2, lng2 = list(map(math.radians, (lat1, lng1, lat2, lng2)))
        
                # calculate haversine
                lat = lat2 - lat1
                lng = lng2 - lng1
                d = math.sin(lat * 0.5) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(lng * 0.5) ** 2
                h = 2 * AVG_EARTH_RADIUS * math.asin(math.sqrt(d))
                if miles:
                    return h * MILES_PER_KILOMETER # in miles
                else:
                    return h # in kilometers 
            Node2Len = dict(list(zip(list(range(len(Point2Node)+len(Segs))), [0]*len(Point2Node)+[haversine(s[1],s[0]) for s in Segs])))

        #Fill the i and j arrays for sparse COO format. 
        [list(map(i.extend, [[Point2Node[l[ind]], Seg2Node[(l[ind], l[ind+1])], Seg2Node[(l[ind], l[ind+1])], Point2Node[l[ind+1]]] for ind in range(len(l)-1)])) for l in Ltup if l !=[]]
        [list(map(j.extend, [[Seg2Node[(l[ind], l[ind+1])], Point2Node[l[ind]], Point2Node[l[ind+1]], Seg2Node[(l[ind], l[ind+1])]] for ind in range(len(l)-1)])) for l in Ltup if l != []]
        AdjDict = sp.coo_matrix((np.ones(len(i)), (i,j))) #The adjacency matrix
        return AdjDict, Seg2Node, Node2Len


    def BuildEdgeIncidence(self, Adj):
        '''Build an oriented edge-incident matrix from an adjacency matrix'''
        UpTriAdj = sp.triu(Adj)
        numEdge = UpTriAdj.nnz
        i = np.array(list(range(numEdge)) + list(range(numEdge)))
        j = np.hstack([UpTriAdj.nonzero()[0], UpTriAdj.nonzero()[1]])
        data = np.array([1.0]*numEdge + [-1.0]*numEdge)
        return sp.coo_matrix((data, (i,j))) #The oriented edge-incidence matrix

    def BigConnComp(self, AdjMat, Point2Node, Seg2Node, Node2Len):
        '''Returns the adjacency matrix of the largest connected component 
        in AdjMat. In other words, we filter for the largest connected 
        component and discard the remaining geometric network'''
        (numComp, ConComp) = sp.csgraph.connected_components(AdjMat)
        numPoints = AdjMat.shape[0]-len(Seg2Node)
        if numComp == 1:
             return (AdjMat, Point2Node, Seg2Node, Node2Len)
        #Find the big connected component and its adjacency matrix
        BigCompInd = np.argmax([sum(ConComp == ind) for ind in range(numComp)])
        BigComp = (ConComp == BigCompInd)
        NewAdjMat = AdjMat[:, BigComp]
        NewAdjMat = NewAdjMat[BigComp, :]
        
        #Remove the removed nodes form Seg2Node and Point2Node
        Node2Seg = dict(list(zip(list(Seg2Node.values()), list(Seg2Node.keys()))))
        [Node2Seg.pop(n) for n in list(Seg2Node.values()) if not BigComp[n]]
        Node2Point = dict(list(zip(list(Point2Node.values()), list(Point2Node.keys()))))
        [Node2Point.pop(n) for n in list(Point2Node.values()) if not BigComp[n]]
           
        #Update the dictionaries 
        Point2Node = dict(list(zip(list(Node2Point.values()), list(range(len(Node2Point))))))
        Seg2Node = dict(list(zip(list(Node2Seg.values()), list(range(len(Node2Point), NewAdjMat.shape[0])))))
        Node2Len = dict(list(zip(list(range(NewAdjMat.shape[0])), [0]*NewAdjMat.shape[0])))
        Node2Len.update(dict(list(zip(list(range(len(Node2Point), NewAdjMat.shape[0])), [np.linalg.norm(np.subtract(i[1],i[0])) for i in list(Seg2Node.keys())]))))
          
        return (NewAdjMat, Point2Node, Seg2Node, Node2Len)

    def SolveProblem(self, lam, solver = "osqp", dual = True, eps = 1e-4):
        '''Solve a previously generated fused density estimation problem for 
        a penalty parameter lam. The default setting is to solve the dual 
        problem to 1e-4 accuracy using the osqp solver.

        Example: 
        fde = FDE(L,P)
        fde.GenerateProblem()
        fde.SolveProblem(.04)  <-- Solves problem for lambda = .04
        fde.SolveProblem(.08)  <-- Easily solve for different penalty'''
        if self.Problem_Declared == False:
            raise ValueError("Problem must be generated before it can be solved")
        if dual == True:
            if solver == "osqp":
                start = timer()
                m = osqp.OSQP()
                m.setup(P = self.D1*sp.diags(1/self.s)*self.D1.T, q = np.zeros(self.D1.shape[0]), A = sp.vstack([self.D2.T, sp.eye(self.D2.shape[0])]), l = np.hstack([-self.u, -lam*np.ones(self.D2.shape[0])]), u = np.hstack([-self.u, lam*np.ones(self.D2.shape[0])]), eps_abs = 1e-4, max_iter = 100000, warm_start = False)
                results = m.solve()
                end = timer()
                #print("OSQP: " + str(end-start))
                if results.info.status_val == 1:
                    self.z = -np.squeeze(self.D1.T*results.x)/self.s
                else:
                    print("Warning! FDE not solved. Increase lambda parameter")
                    self.z = np.zeros(self.D1.shape[1])
    
            if solver == "cvxopt":
                P_coo = (self.D1*sp.diags(1/self.s)*self.D1.T).tocoo()
                P_sparse = cvxopt.spmatrix(list(map(float, P_coo.data)), P_coo.row.tolist(), P_coo.col.tolist(), size = P_coo.shape)
                q = cvxopt.matrix(np.zeros(self.D1.shape[0]))
                A_coo = (self.D2.T).tocoo()
                A_sparse = cvxopt.spmatrix(list(map(float, A_coo.data)), A_coo.row.tolist(), A_coo.col.tolist(), size = A_coo.shape)
                b = cvxopt.matrix(-self.u)
                G_coo = (sp.vstack([sp.eye(self.D2.shape[0]), - sp.eye(self.D2.shape[0])])).tocoo()
                G_sparse = cvxopt.spmatrix(list(map(float, G_coo.data)), G_coo.row.tolist(), G_coo.col.tolist(), size = G_coo.shape)
                h = cvxopt.matrix(np.hstack([lam*np.ones(self.D2.shape[0]), lam*np.ones(self.D2.shape[0])]))
                start = timer()
                results = cvxopt.solvers.qp(P_sparse, q, G_sparse, h, A_sparse, b, options = {'abstol': 1e-4, 'maxiters': 1000})
                end = timer()
                #print("cvxopt: " + str(end-start))
                self.z = -np.squeeze(self.D1.T*results['x'])/self.s
        if dual == False:
            #Possible primals solvers are: cvx.SCS, cvx.ECOS, cvx.CVXOPT, cvx.GUROBI (if installed)
            #start = timer()
            c = cvx.Variable(len(self.s))
            p = cvx.Variable(len(self.u))
            obj = cvx.Minimize(self.u*p+.5*cvx.quad_form(c, cvx.diag(self.s))+lam*cvx.norm(self.D1*c+self.D2*p, 1))
            prob = cvx.Problem(obj)
            if solver == cvx.SCS:
                prob.solve(solver = solver, verbose = True, max_iters = 10000, eps = 1e-5)
            else:
                prob.solve(solver = solver, verbose = True, max_iters = 10000, abstol = 1e-5)
            #end = timer()
            #print(str(solver) + " " + str(end-start))
            self.z = np.squeeze(c.value).tolist()[0]
        self.ScoreVec = np.array([max([self.z[self.D1[j,:].nonzero()[1][0]] for j in self.D2[:,i].nonzero()[0]]) for i in range(self.D2.shape[1])])
        self.Problem_Solved = True
        return end-start

    def score(self, inds):
        if self.Problem_Solved == False:
            raise ValueError("Problem must be solved in order to evaluate score")
 
    # Min u.T * z_2 + norm(D_1 z_1 + D_2 z_2, 1)
    # s.t. S z_1 = -D_1.T *y  <-- y is solved for
    #OR
    # Min u.T * z_2 + 1.T*c
    # c >= D_1 z_1 + D_2 z_2
    # c >= -D_1 z_1 - D_2 z_2
    # S z_1 = -D_1.T*y <-- y is still solved for
    #OR
    #Set z_1 to be the maximum of incident edges.
    #z_1[i] = max(z_2[D2[:,i].nonzero()])
    #z_1[i] = D1 z_

    def fit(self, TrainInds, ValInds, lam = 1):
        if self.Problem_Declared == False:
            raise ValueError("Problem must be declared in order to fit lambda parameter")
        boolInds = np.array([i in ValInds for i in range(len(self.u))])
        self.u = np.copy(self.uTrue)
        self.u[ValInds] = 0 
        self.u = -self.u/sum(self.u)
        self.SolveProblem(lam)
        return self.ScoreVec[ValInds]

    def CrossValidate(self, fold = 20, max_lam =.1):
        '''Lowering max_lam will allow a more precise cross-validation'''
        if self.Problem_Declared == False:
            raise ValueError("Problem must be declared in order to fit lambda parameter")
        self.uTrue = np.copy(self.u)
        Data = list(range(len(self.u)))
        random.shuffle(list(range(len(self.u))))
        PartData = [Data[i:i+len(Data)//fold] for i in range(0, len(Data), len(Data)//fold)]
        Lam = np.linspace(-np.min(self.u)/2, max_lam, 20) #Sensitive
        Results = np.array([[self.fit([i for j in range(fold) if j != k for i in PartData[j]], PartData[k], lam = Lam[l]) for k in range(fold)] for l in range(1, len(Lam))])
        self.u = np.copy(self.uTrue)
        return np.median(Lam[Results.argmax(axis = 0)])

    def Plot(self):
        '''Plot a 2-dimensional fused density estimator.
         If the problem has not been solved, the geometric network and 
        observations will be plotted without the corresponding FDE. 
        Run plt.show() to show results.

        Example:
        fde = FDE(L,P)  <--2D data
        fde.GenerateProblem()
        fde.SolveProblem(.04)
        fde.Plot_2D()'''

        if self.Problem_Declared == False:
            raise ValueError("Problem must be generated in order to plot network")
        ax = plt.axes()
        Segs = np.array(list(self.Seg2s_index.keys()))
        ax = plt.axes()
        ax.set_ylim(np.min(Segs[:,:,1]), np.max(Segs[:,:,1]))
        ax.set_xlim(np.min(Segs[:,:,0]), np.max(Segs[:,:,0]))
        if self.Problem_Solved == True:
            lines = LineCollection(Segs, cmap = "rainbow", array = np.array([self.z[self.Seg2s_index[seg]] for seg in list(self.Seg2s_index.keys())]))
            plt.colorbar(lines, format = '%.2f')
        else:
            lines = LineCollection(Segs, colors = 'b')
        lines.set_clim(vmin = 0)
        ax.add_collection(lines)
        plt.xticks([])
        plt.yticks([])
        plt.scatter([P[0] for P in list(self.Point2Obs.keys()) if P != False and self.Point2Obs[P] != 0] , [P[1] for P in list(self.Point2Obs.keys()) if P != False and self.Point2Obs[P] != 0], s = 50, c = 'k', marker = 'o', zorder = 2)

class UnivarFDE(FDE):
    def __init__(self, xxx_todo_changeme, P):
        (a,b) = xxx_todo_changeme
        if P.min() < a or P.max() > b:
            ValueError("Observations do not lie in interval")
        L = np.unique(np.hstack([[a], P, [b]]))
        L = [[[l] for l in L]]
        FDE.__init__(self, L, P, lat_lon = False)

    def Plot(self):
        '''Plot univariate fused density estimator. Run plt.show()
         to see the results'''
        if self.Problem_Declared == False:
            raise ValueError("Problem must be generated in order to plot network")
        ax = plt.axes()
        if self.Problem_Solved == True:
            #Collect and sort the points
            sortPoints = np.array([[seg[i][0], self.z[self.Seg2s_index[seg]]] for i in [0,1] for seg in list(self.Seg2s_index.keys())])
            sortPoints = sortPoints[sortPoints[:,0].argsort()]
            sortPoints = np.vstack([[sortPoints[0,0], 0], sortPoints, [sortPoints[-1,0], 0]])
            
            horzsegs = np.array([[[x, min(sortPoints[sortPoints[:,0] == x,1])], [x, max(sortPoints[sortPoints[:,0] ==x, 1])]] for x in np.unique(sortPoints[:,0])])
            vertsegs = np.array([[[min(sortPoints[sortPoints[:,1] == y, 0]), y], [max(sortPoints[sortPoints[:,1] == y, 0]), y]] for y in np.unique(sortPoints[:,1])])

            horzlines = LineCollection(horzsegs)
            vertlines = LineCollection(vertsegs)
            ax.add_collection(horzlines)
            ax.add_collection(vertlines)
            ax.set_xlim(np.min(vertsegs[:,:,0]), np.max(vertsegs[:,:,0]))
        P = np.array([p for p in list(self.Point2Obs.keys()) if p != False and self.Point2Obs[p] !=0])
        plt.scatter(P[:], np.zeros(len(P)), s = 50, c = 'r', marker = 'v')
        ax.set_ylim(bottom = 0)

    
