from KFBase import KFVector 
from KFBase import KFMatrix,KFMatrixNull,KFMatrixUnitary 
from kfilter import KFData,KFNode,KFPropagate,KFFilter
from nextphysdat import NEXT,MS,momentum
from math import *
import random

"""
KalmanFilter implementation for NEXT

Trayectory is a straight line (x,y,tx,ty,ene) and the medium is continuous

"""

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
        x,y,tx,ty,ene = x
        p = momentum(ene)
        Qms = MS.QMatrix(p,dz,tx,ty,X0=self.radlen)
        for i in range(4):
            for j in range(4):
                Q[i,j] = Qms[i,j]
        for i in range(4,n): Q[i,i]=1.e-32
        #print ' Q ',Q
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
            
class KFLineFilter(KFFilter):
    """ a KFFilter class for NEXT
    """

    def __init__(self,radlen=0.,ndim=5):
        """ constructor with the xenon presure
        """
        self.hmatrix = KFMatrixNull(2,5)
        self.hmatrix[0,0]=1
        self.hmatrix[1,1]=1
        self.propagator = KFLinePropagate(radlen)
        return

    def sethits(self,hits):
        """ set the hits (KFDatas) with the measurements 
        """
        nodes = map(lambda hit: KFNode(hit,self.hmatrix),hits)
        self.nodes = nodes
        self.status = 'hits'
        return
     
#---- checks -------

def check_kf(radlen=1500.):
    nhits = 4
    
    zdis = 1. # distance z between planes
    z0 = 1. # position of the 1st plane
    tx = 0. # tx slope
    ty = 0. # ty slope
    
    xs = map(lambda i: tx*i*zdis,range(nhits)) # true xs positions
    ys = map(lambda i: ty*i*zdis,range(nhits)) # true ys positions
    zs = map(lambda i: i*zdis+z0,range(nhits)) # true zs positions
    xres = 0.01 # x resolution
    yres = 0.01 # y resolution

    # H matrix 
    V = KFMatrixNull(2) 
    V[0,0]=xres*xres
    V[1,1]=yres*yres

    # hits
    hits = map(lambda x,y,z: KFData(KFVector([x,y]),V,z),xs,ys,zs)
    print ">>> Preparation "
    print ' empty hits ',hits

    # create NEXT Kalman Filter
    kf = KFLineFilter(radlen)
    kf.sethits(hits)
    x = KFVector([0.,0.,0.,0.,2.5])
    C = KFMatrixNull(5,5)
    state0 = KFData(x,C,0.)

    print ">>> Generation "
    hits,vstates = kf.generate(state0)  # generate
    for i in range(nhits):
        print ' state ',vstates[i].vec
        print ' hit ',hits[i]
    
    print '>>>> Fit '
    kf.sethits(hits)
    kf.fit(state0)
    for node in kf.nodes:
        print 'NODE ',node
        print 'chi2 ',node.getchi2('smooth')

#--- main ---    
   
if __name__ == '__main__':

    check_kf()
