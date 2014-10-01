from KFBase import KFVector 
from KFBase import KFMatrix,KFMatrixNull,KFMatrixUnitary 
from kfilter import KFData,KFNode,KFPropagate,KFFilter
from math import *
import random

"""
KalmanFilter implementation for NEXT

Trayectory is a straight line (x,y,tx,ty,ene) and the medium is continuous

"""

def momentum(ene, mass = 0.511):
    """ return the momentum from a given energy 
    """
    p = sqrt((ene+mass)**2-mass**2)
    return p

def thetams(p, dis, mass= 0.511):
    """ return the theta MS angle for a particle with mass, momentum (p) that travel a distance (dis)
    """
    ene = sqrt(mass**2+p**2)
    beta = p/ene
    tms = (13.6)/(p*1.*beta)*sqrt(dis*1.)*(1+0.038*log(dis*1.))
    return tms
    

class KFLinePropagate(KFPropagate):
    """ it assumes and state (x,y,tx,ty,ene) 
    """

    def __init__(self,radlen=0.):
        """ constructor of a LinePropoagator with a redlen
        """
        self.radlen = radlen
        return

    def QMatrix(self,x,dz):
        """ Returns the Q matrix for a deltaz
        """ 
        n = x.Length()
        Q = KFMatrixNull(n)
        if (self.radlen == 0.): return Q
        x,y,tx,ty,ene = svec[0],svev[1],svec[2],svec[3],svec[4]
        tx2 = tx*tx
        ty2 = ty*ty
        tt2 = 1.+tx2+ty2
        ds = sqrt(tt2)*abs(dz)
        dsig = 1.
        if (dz<0): dsig = -1.
        rad = ds/self.radlen
        pp = momentum(ene)
        tms = thetams(p,rad)
        print ' rad, pp, thetams ',rad,pp,tms
        norm2 = tt2*tms*tms
        covtxtx = norm2*(1+tx2) 
        covtyty = norm2*(1+ty2) 
        covtxty = norm2*tx*ty
        Q[0,0] = ds2*covtxtx/3.
        print ' ds covtx covty covtxty norm ',ds,covtxtx,covtyty,covtxty,norm
        Q[1,0] = ds2*covtxty/3.
        Q[2,0] = covtxtx*ds*dsign
        Q[3,0] = covtxty*ds*dsign
        Q[1,1] = ds2*covtyty/3.
        Q[2,1] = covtxty*ds*dsign
        Q[3,2] = covtyty*ds*dsign
        Q[2,2] = covtxtx
        Q[2,3] = covtxty
        Q[3,3] = covtyty
        for i in range(n):
            for j in range(i+1,n):
                Q[i,j] = Q[j,i]        
        print ' QMatrix ',Q
        return Q

    def FMatrix(self,x,dz):
        """ return the FMatrix
        """
        n = x.Length()
        F = KFMatrixUnitary(n)
        F[0,2]=dz
        F[1,3]=dz
        #print ' FMatrix ',F
        return F

    def propagate(self,state,zrun):
        """ propagate this state to a zrun position
        """
        #print ' propagate state ',state
        #print ' propagate at ',zrun
        x = KFVector(state.vec)
        C = KFMatrix(state.cov)
        deltaz = zrun-state.zrun
        F = self.FMatrix(x,deltaz)
        FT = F.Transpose()
        #print ' F ',F
        #print ' FT ',FT
        Q = self.QMatrix(x,deltaz)
        #print ' Q ',Q
        xp = F*x
        #print ' xp ',xp
        Cp = F*C*FT+Q
        #print ' Cp ',Cp
        pstate = KFData(xp,Cp,zrun)
        #print ' propagated state ',pstate
        #print " F ",F
        #print " Q ",Q
        return pstate,F,Q
    
class KFNextFilter(KFFilter):
    """ a KFFilter class for NEXT
    """

    def __init__(self,presure=0.,ndim=5):
        """ constructor with the xenon presure
        """
        self.hmatrix = KFMatrixNull(2,5)
        self.hmatrix[0,0]=1
        self.hmatrix[1,1]=1
        radlen = presure
        self.propagator = KFLinePropagate(radlen)
        return

    def sethits(self,hits):
        """ set the hits (KFDatas) with the measurements 
        """
        nodes = map(lambda hit: KFNode(hit,self.hmatrix),hits)
        self.nodes = nodes
        self.status = 'hits'
        return
        
if __name__ == '__main__':

    # Kalman Filter with a Straight Line  

    nhits = 4 # number of hits
    zdis = 1. # distance z between planes
    z0 = 1. # position of the 1st plane
    tx = 1. # tx slope
    ty = 0. # ty slope
    xs = map(lambda i: tx*i*zdis,range(nhits)) # true xs positions
    ys = map(lambda i: ty*i*zdis,range(nhits)) # true ys positions
    zs = map(lambda i: i*zdis+z0,range(nhits)) # true zs positions
    xres = 0.01 # x resolution
    yres = 0.01 # y resolution
    
    xss = map(lambda x: x+random.gauss(0.,xres),xs) # measured xs positions
    yss = map(lambda x: x+random.gauss(0.,yres),ys) # measured xs positions
    print ' xss ',xss
    print ' yss ',yss
    print ' zs ',zs

    # H matrix 
    V = KFMatrixNull(2) 
    V[0,0]=xres*xres
    V[1,1]=yres*yres
    print ' V ',V

    # hits
    hits = map(lambda x,y,z: KFData(KFVector([x,y]),V,z),xss,yss,zs)
    print ' hits ',hits
    
    # create NEXT Kalman Filter
    kf = KFNextFilter()
    kf.sethits(hits)
    x = KFVector([0.,0.,0.,0.,2.5])
    C = KFMatrixUnitary(5)
    state0 = KFData(x,C,0.)
    kf.filter(state0)  # filter
    kf.smoother()

    print '>>>> Results '
    for node in kf.nodes:
        print 'NODE ',node
        print 'chi2 ',node.getchi2('smooth')

    
