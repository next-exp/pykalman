from KFBase import KFVector 
from KFBase import KFMatrix,KFMatrixNull,KFMatrixUnitary,Random 
from kfilter import KFData,KFNode,KFPropagate,KFFilter
from nextkfnoiser import NEXT,MSNoiser,ELoss
from math import *
import random

DEBUG = False

"""
KalmanFilter implementation for NEXT

Trayectory is a straight line (x,y,tx,ty,ene) and the medium is continuous

"""

class KFZLinePropagate(KFPropagate):
    """ StraightLine propagation of a ZState
    """

    def __init__(self,noiser=None,eloss=None):
        """ constructor of a ZLine propoagator for a (x,y,tx,ty,ene) statte with
        para the forward(+1)/backward(-1) sense along z
        It uses a noiser for the MS and energy loss for ene parameters
        """
        self.noiser=noiser
        self.eloss=eloss 
        return

    def validstep(self,state,zrun):
        dz = zrun-state.zrun
        ok = True
        #ok = (state.par*dz>0.)
        if (not ok): 
            print "KFZLinePropagate.velid step Warning no valid dz ",dz,' state z-forward/backward ',state.par
            return ok
        ene = state.vec[-1]
        if (self.noiser):
            ok = self.noiser.validstep(ene,zrun)
            if (not ok):
                print "KFZLinePropagate.velid step Warning no valid noiser "
        if (self.eloss):
            ok = self.eloss.validstep(ene,zrun)
            if (not ok):
                print "KFZLinePropagate.velid step Warning no valid noiser "
        return ok
        
    def QMatrix(self,x,dz):
        """ Returns the Q matrix for a deltaz
        """ 
        n = x.Length()
        Q = KFMatrixNull(n)
        if (not self.noiser): return Q
        x,y,tx,ty,ene = x
        p = ene
        Qms = self.noiser.QMatrix(p,dz,tx,ty)
        for i in range(4):
            for j in range(4):
                Q[i,j] = Qms[i,j]
        for i in range(4,n): Q[i,i]=1.e-32
        if (DEBUG): print 'KFZLinePropagate.QMatrix ',Q
        return Q

    def FMatrix(self,x,dz):
        """ return the FMatrix
        """
        n = x.Length()
        F = KFMatrixUnitary(n)
        F[0,2]=dz
        F[1,3]=dz
        if (self.eloss):
            tx,ty,ene0=x[2],x[3],x[4]
            icost = sqrt(1.+tx*tx+ty*ty)
            ds = dz*icost
            dene,ds = self.eloss.deltaE(ene0,ds)
            rat = (ene0-dene)/ene0
            assert rat >= 0,'KFZLinePropoage.FMatrix no valid eloss %f'%rat
            F[4,4] = rat
        if (DEBUG): 
            print 'KFNextPropagate.FMatrix x ',x,' dz ',dz,' F ',F
        return F

    def propagate(self,state,zrun):
        """ propagate this state to a zrun position
        """
        #print ' propagate state ',state
        #print ' propagate at ',zrun
        ok = self.validstep(state,zrun)
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
        pstate = KFData(xp,Cp,zrun,state.pars)
        #print ' propagated state ',pstate
        #print " F ",F
        #print " Q ",Q
        if (DEBUG): print 'KFZLinePropagate.propagate state,F,Q ',pstate,F,Q
        return ok,pstate,F,Q

class KFZLineFilter(KFFilter):
    """ a KFFilter class for NEXT
    """

    def __init__(self,noiser=None,eloss=None):
        """ constructor with the xenon presure
        """        
        self.hmatrix = KFMatrixNull(2,5)
        self.hmatrix[0,0]=1
        self.hmatrix[1,1]=1
        self.propagator = KFZLinePropagate(noiser,eloss)
        return

    def sethits(self,hits):
        """ set the hits (KFDatas) with the measurements 
        """
        nodes = map(lambda hit: KFNode(hit,self.hmatrix),hits)
        self.nodes = nodes
        self.status = 'hits'
        return

    def setnodes(self,nodes):
        self.nodes = nodes
        self.state = 'hits'
        return
     
class KFNextGenerator:

    def __init__(self,deltae=0.05,deltax=0.1):
        self.deltae = deltae
        self.deltax = deltax
        self.next = NEXT()
        self.msnoiser = self.next.msnoiser
        self.eloss = self.next.eloss
        self.propagator = KFZLinePropagate(self.msnoiser,self.eloss)
        return

    def generate(self,state0):
        states = []
        state = list(state0)
        ee = state[-1]
        while (self.eloss.validstep(ee,self.deltax)):
            states.append(state)
            x0,y0,z0,ux0,uy0,uz0,ee0 = state
            xv0 = KFVector([x0,y0,z0])
            uv0 = KFVector([ux0,uy0,uz0])
            uv0.Unit()
            #print ' x0 ',xv0
            #print ' u0 ',uv0
            #print ' ee0 ',ee0
            de,ds = self.eloss.deltax(ee0,self.deltae)
            #print ' de, ds ',de,ds
            xvf,uvf = self.msnoiser.XUrandom(ee0,ds,xv0,uv0)
            x,y,z=xvf; ux,uy,uz=uvf
            ee = ee0-de
            state = [x,y,z,ux,uy,uz,ee]
            if (DEBUG): 
                print 'KFNextGenerator.generate new state ',state
        return states

    @staticmethod
    def kfnodes(states,zs,xres):
        H = KFMatrixNull(2,5); H[0,0]=1.; H[1,1]=1.
        zsts = KFNextGenerator.sample(states,zs)
        hits = KFNextGenerator.hits(zsts,xres)
        kfts = KFNextGenerator.kfstates(zsts)
        def knode(hit,st):
            nod = KFNode(hit,H)
            nod.setstate('true',st)
            return nod
        ks = map(knode,hits,kfts)
        if (DEBUG):
            print "KFNextGenerator.kfnodes kfnodes ",ks
        return ks

    @staticmethod
    def sample(states,zs):
        zstates = []
        nn = len(states)
        def getzi(st0,st1):
            x0,y0,z0,ux0,uy0,uz0,ee0=st0
            x1,y1,z1,ux1,uy1,uz1,ee1=st1
            for zi in zs:
                ok = False
                if (z0<=zi and zi<z1): ok =True
                if (z0>=zi and zi>z1): ok = True
                if (ok):
                    #print ' zi! ',z0,zi,z1
                    return zi
            return None
        def zstate(st0,st1,z):
            x0,y0,z0,ux0,uy0,uz0,ee0=st0
            x1,y1,z1,ux1,uy1,uz1,ee1=st1
            x = x0+(ux0/uz0)*(z-z0)
            y = y0+(uy0/uz0)*(z-z0)
            ee = ee0-(ee0-ee1)*(z-z0)/(z1-z0)
            stz = (x,y,z,ux0,uy0,uz0,ee)
            #print " st0 ",st0
            #print " st1 ",st1
            #print " stz ",stz
            return stz
        for ii in range(1,nn-1):
            st0,st1 = states[ii],states[ii+1]
            zi = getzi(states[ii],states[ii+1])
            if (not zi): continue
            zst =  zstate(st0,st1,zi)
            zstates.append(zst)
        if (DEBUG): 
            print "KFNextGenerator.sample zstates ",zstates
        return zstates
    
    @staticmethod
    def hits(states,xres=0.01):
        zhits = []
        V = KFMatrix([[xres*xres,0.],[0.,xres*xres]])
        for state in states:
            x0,y0,z0,ux,uy,uz,ee = state            
            m0 = KFVector([x0,y0])
            sm = Random.cov(V)
            mm = m0+sm
            zhit = KFData(mm,V,z0)
            zhits.append(zhit)
        if (DEBUG): 
            print "KFNextGenerator.hits hits ",zhits
        return zhits

    @staticmethod
    def kfstates(states):
        sts = []
        C0 = KFMatrixNull(5,5)
        for state in states:
            x,y,z,ux,uy,uz,ee = state
            st = KFData(KFVector([x,y,ux/uz,uy/uz,ee]),C0,z,pars={'uz':uz})
            sts.append(st)
        if (DEBUG):
            print "KFNextGenerator.kfstates ",sts
        return sts

#---- checks -------

def ck_kfilter(radlen=1500.):
    
    # Preparation
    #--------------

    # V matrix 
    xres = 0.01 # x resolution
    yres = 0.01 # y resolution
    V = KFMatrixNull(2) 
    V[0,0]=xres*xres
    V[1,1]=yres*yres

    # hits
    nhits = 10 # number of hits
    zdis = 2. # distance z between planes
    z0 = 1. # position of the 1st plane
    zs = map(lambda i: z0+zdis*i,range(nhits))
    hits = map(lambda z: KFData(KFVector([0.,0.]),V,z),zs)
    print ">>> Preparation "
    print ' empty hits ',hits

    # create NEXT Kalman Filte
    msnoiser = MSNoiser(radlen)
    kf = KFZLineFilter(noiser=msnoiser)
    kf.sethits(hits)

    print ">>> Generation "
    x0 = KFVector([0.,0.,0.,0.,2.5])
    C0 = KFMatrixNull(5,5)
    state0 = KFData(x0,C0,0.)
    knodes = kf.generate(state0)  # generate
    for i in range(len(kf)):
        print 'knode ',i,knodes[i]
        
    print '>>>> Fit '
    x0 = KFVector([0.,0.,0.,0.,2.5])
    C0 = KFMatrixUnitary(5)
    state0 = KFData(x0,C0,0.)
    kf.setnodes(knodes)
    ok,fchi2,schi2 = kf.fit(state0)
    print "ok, fchi2, schi2 ",ok,fchi2,schi2
    if (not ok):
        print "Not possible to fit! "
        return
    for node in kf.nodes:
        print 'node ',node

def ck_generator():

    from alex.rootsvc import ROOTSvc
    root = ROOTSvc('root','temp.root')
    root.open()

    deltae=0.01
    nhits = int(2.5/deltae)

    root.h2d('xi',nhits,0,1.*nhits,100,-15.,15.)
    root.h2d('yi',nhits,0,1.*nhits,100,-15.,15.)
    root.h2d('zi',nhits,0,1.*nhits,100,-5.,20.)
    root.h2d('ei',nhits,0,1.*nhits,100,0.,3.)
    
    root.h2d('xz',100,-5.,20,100,-15.,15.)
    root.h2d('yz',100,-5.,20,100,-15.,15.)

    root.h2d('kxz',100,-5.,20,100,-15.,15.)
    root.h2d('kyz',100,-5.,20,100,-15.,15.)

    for k in range(10):

        pp = '_evt'+str(k)

        root.h2d('xi'+pp,nhits,0,1.*nhits,100,-15.,15.)
        root.h2d('yi'+pp,nhits,0,1.*nhits,100,-15.,15.)
        root.h2d('zi'+pp,nhits,0,1.*nhits,100,-5.,20.)
        root.h2d('ei'+pp,nhits,0,1.*nhits,100,0.,3.)
        
        root.h2d('xz'+pp,100,-5.,20,100,-15.,15.)
        root.h2d('yz'+pp,100,-5.,20,100,-15.,15.)

        root.h2d('kxz'+pp,100,-5.,20,100,-15.,15.)
        root.h2d('kyz'+pp,100,-5.,20,100,-15.,15.)
    
        kfgen = KFNextGenerator(deltae)

        state0 = [0.,0.,0.,0.,0.,1.,2.5]
        states = kfgen.generate(state0)
        print " states ",states

        dz = 0.5
        zs = map(lambda i: (i+1)*dz,range(20))
        #zstates = KFNextGenerator.sample(states,zs)
        #print " zstates ",zstates

        #hits = KFNextGenerator.hits(zstates,0.01)
        #print " zstates ",hits

        nods = KFNextGenerator.kfnodes(states,zs,0.01)
        print " znodes ",nods

        i = 0
        for state in states:
            i+=1
            x,y,z,ux,uy,uz,ee = state
            root.fill('xi'+pp,i,x)
            root.fill('yi'+pp,i,y)
            root.fill('zi'+pp,i,z)
            root.fill('ei'+pp,i,ee)
            root.fill('xz'+pp,z,x)
            root.fill('yz'+pp,z,y)

            root.fill('xi',i,x)
            root.fill('yi',i,y)
            root.fill('zi',i,z)
            root.fill('ei',i,ee)
            root.fill('xz',z,x)
            root.fill('yz',z,y)

        for nod in nods:
            st = nod.getstate('true')
            x,y,tx,ty,ee = st.vec; z = st.zrun
            root.fill('kxz'+pp,z,x)
            root.fill('kyz'+pp,z,y)

            root.fill('kxz',z,x)
            root.fill('kyz',z,y)

    root.close()
#--- main ---    
   
if __name__ == '__main__':

    #ck_kfilter(1500.)
    ck_generator()
