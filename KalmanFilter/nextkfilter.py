from KFBase import KFVector 
from KFBase import KFMatrix,KFMatrixNull,KFMatrixUnitary,Random 
from kfilter import KFData,KFNode,KFPropagate,KFFilter
from nextphysdat import NEXT,ELoss,MS
from math import *
import random

DEBUG = False

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
        self.ms = MS(radlen)
        return

    def isvalid(self,state,zrun):
        xx = state.vec
        x0,y0,tx,ty,ee = xx
        z0 = state.zrun
        icost = sqrt(1+tx*tx+ty*ty)
        ds = (zrun-z0)*icost
        if (abs(ds) == 0.): return False
        ok = self.ms.isvalid(ee,abs(ds))
        if (not ok): 
            print "Warning! propagate no valid due to ms"
            return ok
        return ok

    def QMatrix(self,x,dz):
        """ Returns the Q matrix for a deltaz
        """ 
        n = x.Length()
        Q = KFMatrixNull(n)
        if (self.radlen == 0.): return Q
        x,y,tx,ty,ene = x
        p = ene
        Qms = self.ms.QMatrix(p,dz,tx,ty)
        for i in range(4):
            for j in range(4):
                Q[i,j] = Qms[i,j]
        for i in range(4,n): Q[i,i]=1.e-32
        if (DEBUG): print 'KFLinePropagate.QMatrix ',Q
        return Q

    def FMatrix(self,x,dz):
        """ return the FMatrix
        """
        n = x.Length()
        F = KFMatrixUnitary(n)
        F[0,2]=dz
        F[1,3]=dz
        if (DEBUG): print 'KFLinePropagate.FMatrix ',F
        return F

    def propagate(self,state,zrun):
        """ propagate this state to a zrun position
        """
        #print ' propagate state ',state
        #print ' propagate at ',zrun
        ok = self.isvalid(state,zrun)
        if (not ok): 
            print "Warning! not valid propagate at ",zrun," state ",state.vec
            return ok,None,None,None
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
        if (DEBUG): print 'KFLinePropagate.propagate state,F,Q ',pstate,F,Q
        return ok,pstate,F,Q
            
class KFNextPropagate(KFLinePropagate):

    def __init__(self,next=None):
        self.next = next
        if (not self.next): self.next = NEXT()
        self.ms = self.next.ms
        self.eloss = self.next.eloss
        self.radlen = self.ms.X0
        return

    def isvalid(self,state,zrun):
        ok = KFLinePropagate.isvalid(self,state,zrun)
        if (not ok): return ok
        xx = state.vec
        x0,y0,tx,ty,ee = xx
        z0 = state.zrun
        icost = sqrt(1+tx*tx+ty*ty)
        ds = abs(zrun-z0)*icost
        ok = self.eloss.isvalid(ee,ds)
        if (not ok): 
            print "Warning! propagate no valid due to eloss"
            return ok
        return True
        
    def FMatrix(self,x,zrun):
        F = KFLinePropagate.FMatrix(self,x,zrun)
        tx,ty,ene0=x[2],x[3],x[4]
        if (ene0<=0.): self.success = False
        icost = sqrt(1+tx*tx+ty*ty)
        ds = zrun*icost
        dene,ds = self.eloss.deltaE(ene0,ds)
        F[4,4] = (ene0-dene)/ene0
        if (DEBUG): 
            print 'KFNextPropagate.FMatrix ene0,dene,ds,F',ene0,dene,ds,F
        return F

class KFLineFilter(KFFilter):
    """ a KFFilter class for NEXT
    """

    def __init__(self,eloss=False,radlen=-1.):
        """ constructor with the xenon presure
        """
        self.hmatrix = KFMatrixNull(2,5)
        self.hmatrix[0,0]=1
        self.hmatrix[1,1]=1
        next = NEXT()
        if (not eloss):
            if (radlen<0.): radlen = next.X0
            print 'setting propagator with radlen ',radlen 
            self.propagator = KFLinePropagate(radlen)
        else:
            print 'setting propagator next radlen ',next.ms.X0
            self.propagator = KFNextPropagate(next)
        return

    def sethits(self,hits):
        """ set the hits (KFDatas) with the measurements 
        """
        nodes = map(lambda hit: KFNode(hit,self.hmatrix),hits)
        self.nodes = nodes
        self.status = 'hits'
        return

    #def user_filter(self,node):
        #fstate = node.getstate('filter')
        #tstate = node.getstate('true')
        #fstate.vec[-1] = tstate.vec[-1]
    #    return
     
class KFNextGenerator:

    def __init__(self,deltae=0.05):
        self.mass = 0.511 
        self.deltae = deltae
        self.next = NEXT()
        self.ms = MS(self.next.X0)
        self.eloss = ELoss(self.next.rho)
        # self.propagator = KFLinePropagator(radlen=self.X0)
        self.propagator = KFNextPropagate()
        return

    def generate(self,state0,H,V):
        knodes = []
        state = state0.copy()
        C0 = KFMatrixNull(5,5)
        zprev = state.zrun
        while (state.vec[-1]>self.deltae):
            xx = state.vec
            x0,y0,tx0,ty0,ee0 = xx
            z0 = state.zrun
            m0 = H*xx
            #print ' m0 ',m0
            sm = Random.cov(V)
            m = m0+sm
            #print ' hit ',m
            hit = KFData(m,V,z0)
            knode = KFNode(hit,H)
            knode.setstate('true',state)
            knodes.append(knode)
            xv0 = KFVector([x0,y0,z0]) 
            uz = 1.
            if (zprev>z0): uz=-1.
            zprev = z0
            xv0 = KFVector([x0,y0,z0])
            uv0 = KFVector([tx0,ty0,uz])
            uv0.Unit()
            #print ' x0 ',xv0
            #print ' u0 ',uv0
            de,ds = self.eloss.deltax(ee0)
            #print ' de, ds ',de,ds
            p = ee0
            xvf,uvf = self.ms.XUrandom(p,ds,xv0,uv0)
            #print ' xf, vf ',xvf,uvf
            xf,yf,zf = xvf; uxf,uyf,uzf = uvf; txf = uxf/uzf; tyf = uyf/uzf
            xx = KFVector([xf,yf,txf,tyf,ee0-de])
            #print ' >> xx ',xx,' z >> ',zf
            state = KFData(xx,C0,zf)
            #print ' state ',state
        return knodes

class KFNextLenghtGenerator:

    def __init__(self,deltae=0.05):
        self.mass = 0.511 
        self.deltae = deltae
        self.next = NEXT()
        self.ms = MS(self.next.X0)
        self.eloss = ELoss(self.next.rho)
        # self.propagator = KFLinePropagator(radlen=self.X0)
        self.propagator = KFNextPropagate()
        return

    def generate(self,state0,H,V):
        knodes = []
        state = state0.copy()
        C0 = KFMatrixNull(5,5)
        zprev = state.zrun
        while (state.vec[-1]>self.deltae):
            # store states
            xx = state.vec
            x0,y0,z0,tx0,ty0,ee0 = xx
            z0 = state.zrun
            m0 = H*xx
            sm = Random.cov(V)
            m = m0+sm
            hit = KFData(m,V,z0)
            knode = KFNode(hit,H)
            knode.setstate('true',state)
            knodes.append(knode)
            # generate new state
            de,ds = self.eloss.deltax(ee0)
            icost = sqrt(1+tx*tx+ty*ty)
            dz = ds/icost
            ok,state,F,Q = self.propagator(state,dz)
            if (not ok): break
        if (DEBUG):
            print 'KFNextLengthGenerator.propagate ',knodes
        return knodes

#---- checks -------

def ck_kfilter(radlen=1500.):
    nhits = 10
    
    zdis = 2. # distance z between planes
    z0 = 1. # position of the 1st plane
    tx = 1. # tx slope
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
    for i in range(len(hits)):
        print ' state ',vstates[i].vec
        print ' hit ',hits[i]
    
    print '>>>> Fit '
    kf.sethits(hits)
    ok,fchi2,schi2 = kf.fit(state0)
    print "ok, fchi2, schi2 ",ok,fchi2,schi2
    if (not ok):
        print "Not possible to fit! "
        return
    for node in kf.nodes:
        print 'NODE ',node
        print 'chi2 ',node.getchi2('smooth')

def ck_generator():

    E0 = 2.5; z0 = 0.
    x0 = KFVector([0.,0.,0.,0.,E0])
    C0 = KFMatrixNull(5,5)
    V = KFMatrixNull(3,3)
    xres = 0.01; eres=0.01
    V[0,0]=xres*xres; V[1,1]=xres*xres; V[2,2]=eres*eres
    H = KFMatrixNull(3,5)
    H[0,0]=1.;H[1,1]=1.;H[2,4]=1.;
    state0 = KFData(x0,C0,z0)

    kfgen = KFNextGenerator()

    hits,states = kfgen.generate(state0,H,V)

    print "Hits ",hits
    print "States ",states


#--- main ---    
   
if __name__ == '__main__':

    ck_kfilter(1500.)
    #hits,states = ck_generator()
