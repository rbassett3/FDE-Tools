import xml.etree.ElementTree   
import numpy as np
import io
import os
import time

def LoadMapXML(fpath):
    '''Load an openmap xml file.

    LoadMapXML takes a filepath as input and returns a list of numpy arrays.
    Each array corresponds to a highway in the map, and is of size n x 2, 
    where n in the number of nodes in the way. The first column is longitude
    (E-W) of the point and the 2nd is latitude (N-S).'''
    
    f = xml.etree.ElementTree.parse(fpath).getroot()
    W = f.findall('.//way/tag[@k="highway"]/..')                                        
    #xpath syntax is used for following commands
    L = [np.array([map(f.find(''.join(['.//node[@id="', n.get('ref'), '"]'])).get, ["lon", "lat"]) for n in w.findall('nd')]) for w in W]
    L = [l.astype(float) for l in L]
    return L

def GetEateries(fpath):
    '''Get the coordinates of the eateries in a map.

    GetEateries takes a filepath as input and returns a numpy array. Each row corresponds to the coordinate of an eatery. The first column is longitude (E-W) and the second in latitude (N-S).'''

    f = xml.etree.ElementTree.parse(fpath).getroot()
    R = f.findall('.//node/tag[@v="restaurant"]/..')
    R.extend(f.findall('.//node/tag[@v="cafe"]/..'))
    R.extend(f.findall('.//node/tag[@v="fast_food"]/..'))

    return np.array([[float(r.get('lon')), float(r.get('lat'))] for r in R])

def FilterData(D, (lon_lb, lon_ub), (lat_lb, lat_ub)):
    '''Filter the data D, an n x 2 matrix by longitudinal bounds [lon_lb, lon_ub]. In other words, 
    remove data which does not have it longitude greater than lon_lb and less and lon_ub. Do the same 
    with the latitude with respect to lat_lb and lat_ub'''
    #Purge in the following lines
    newInd = np.nonzero(np.logical_and(np.logical_and(D[:,0]>=lon_lb, D[:,0]<=lon_ub), np.logical_and(D[:,1]>=lat_lb, D[:,1]<=lat_ub)))[0]    
    return D[newInd,:]

def DiscretizePaths(L, numpoints = 20): 
    '''Take each segment l in the geometric network L and add numpoints evenly
    spaced points between the points determining l. This allows a more accurate
    projection of observations onto points of the network, but is very 
    expensive.'''
    NewL =[]
    for l in L:
        Newl = l[0,:]
        for ind in range(l.shape[0]-1):
            Newl = np.vstack([Newl, np.array([lam*l[ind,:]+(1-lam)*(l[ind+1,:]) for lam in np.linspace(0,1,numpoints)[1:]])])
        NewL.append(Newl)
    return NewL


