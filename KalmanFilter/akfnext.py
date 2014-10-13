from alex.alex import IAlg
from alex.alex import Alex
from alex.rootsvc import ROOTSvc
import random
from KFBase import KFVector, KFMatrix, KFMatrixNull, KFMatrixUnitary
from kfnext import NEXT, nextgenerator, nextfilter, V0, H0
from kfnext import simplegenerator, simplefilter
from kfgenerator import zsample
from kffilter import randomnode, KFData
from kfzline import zstate
from math import *

"""

Alex example to check KFNextFilter

"""

class Gen(IAlg):

    def define(self):
        self.pgas = 10.
        self.radlen = 1500.
        self.deltae = 0.01 # 10 keV
        self.nhits = 20
        self.dz = 0.5
        self.xres = 0.1
        self.E0 = 2.5
        return

    def execute(self):
        
        next = NEXT(self.pgas)
        kfgen = nextgenerator(next,deltae=self.deltae)
        #kfgen = simplegenerator(self.radlen,deltae=self.deltae)
        state0 = [0.,0.,0.,0.,0.,1.,self.E0]
        nstates = int(self.E0/self.deltae)
        states = kfgen.generate(state0)
        print " states ",len(states)        

        zs = map(lambda i: (i+1)*self.dz,range(self.nhits))        

        V = V0(self.xres)
        zsts = zsample(states,zs)
        print " zstates ",len(zsts)        
        z0 = zsts[0][2]-1.
        for ii,zst in enumerate(zsts):
            if (zst[2] <= z0): break
            z0 = zst[2]
        zsts = zsts[:ii]
        print " zstates forward ",len(zsts)        
        print " state[0] ",zsts[0]
        print " state[-1] ",zsts[-1]

        zzsts = map(zstate,zsts)
        V = V0(self.xres)
        nodes = map(lambda st: randomnode(st,H0,V),zzsts)
        print " nodes ",len(nodes)
        
        ok = False
        if (len(states)>0):
            self.evt['gen/states'] = states
        if (len(nodes)>0.):
            self.evt['gen/nodes'] = nodes
            ok = True

        return ok

class KFGen(IAlg):
    """ Algorithm to generate SL tracks, fit them with KFNextFilter and fill histos
    """

    def define(self):
        self.pgas = 10.
        self.radlen = 1500.
        self.nhits = 20
        self.dz = 0.5
        self.xres = 0.1
        self.E0 = 2.5
        return

    def execute(self):

        # Prepare Generation
        #----------------
        # true zs positions
        zs = map(lambda i: (i+1)*self.dz,range(self.nhits)) 

        # V matrix 
        V = V0(self.xres) 

        # Generate Track
        #----------------
        #empty this
        nullhits = map(lambda z: KFData(KFVector([0,0]),V,z),zs)
        print " null hits ",len(nullhits)
        
        # create NEXT Kalman Filter
        next = NEXT(self.pgas)
        kf = nextfilter(next)
        #kf = simplefilter(self.radlen)
        kf.clear()
        map(lambda hit: kf.addnode(hit,H0),nullhits)
        print " null nodes ",len(kf)

        state0 = zstate([0.,0.,0.,0.,0.,1.,self.E0])
        # generate hits and states into nodes
        knodes = kf.generate(state0)
        print " generated nodes ",len(knodes)
        print " kf nodes ",len(kf)

        ok = (len(knodes)>0)
        if (ok):
            self.evt['gen/nodes']=knodes
        return ok


class KFFit(IAlg):
    """ Algorithm to generate SL tracks, fit them with KFNextFilter and fill histos
    """

    def define(self):
        self.pgas = 10.
        self.radlen = 1500.
        self.E0 = 2.5
        return

    def execute(self):

        nodes = self.evt['gen/nodes']
        if (not nodes): return False

        next = NEXT(self.pgas)
        kf = nextfilter(next)
        #kf = simplefilter(self.radlen)
        kf.clear()
        kf.setnodes(nodes)
        print " fit nodes ",len(nodes) 

        # fit!
        x0 = KFVector([0.,0.,0.,0.,self.E0])
        C1 = KFMatrixUnitary(5)
        state0 = KFData(x0,C1,0.,pars={'uz':1.})
        ok,fchi,schi = kf.fit(state0)  
        if (not ok):
            print " Fit Failed!"
            return False
        print ' Fit chi2 ',fchi,schi

        # set the kalmanfilter
        self.evt['kf'] = kf
        return True

#--------- histogramming!!

class GenHistos(IAlg):

    def define(self):
        self.nevts = 0
        self.nstates = 250
        return
    
    def initialize(self):
        
        self.ievt = -1
        nhits = self.nstates

        self.root.h1d('gh_x',100,-15.,15.)
        self.root.h1d('gh_y',100,-15.,15.)
        self.root.h1d('gh_z',100,-10.,20.)
        self.root.h1d('gh_e',100,0.,3.)

        self.root.h2d('gh_xi',nhits,0,1.*nhits,100,-15.,15.)
        self.root.h2d('gh_yi',nhits,0,1.*nhits,100,-15.,15.)
        self.root.h2d('gh_zi',nhits,0,1.*nhits,100,-5.,20.)
        self.root.h2d('gh_ei',nhits,0,1.*nhits,100,0.,3.)

        self.root.h1d('gh_dx',100,-0.5,0.5)
        self.root.h1d('gh_dy',100,-0.5,0.5)
        self.root.h1d('gh_dz',100,-0.5,0.5)
        self.root.h1d('gh_dth',100,-pi,pi)
        self.root.h1d('gh_dph',100,-pi,pi)
        self.root.h1d('gh_dee',100,-0.2,0.2)

        self.root.h2d('gh_dxi',nhits,0.,1.*nhits,100,-0.5,0.5)
        self.root.h2d('gh_dyi',nhits,0.,1.*nhits,100,-0.5,0.5)
        self.root.h2d('gh_dzi',nhits,0.,1.*nhits,100,-0.5,0.5)
        self.root.h2d('gh_dthi',nhits,0.,1.*nhits,100,-pi,pi)
        self.root.h2d('gh_dphi',nhits,0.,1.*nhits,100,-pi,pi)
        self.root.h2d('gh_deei',nhits,0.,1.*nhits,100,-0.2,0.2)
        
        self.root.h2d('gh_xz',100,-5.,20,100,-15.,15.)
        self.root.h2d('gh_yz',100,-5.,20,100,-15.,15.)
        self.root.h2d('gh_ez',100,-5.,20,100,0.,3.)

    def execute(self):

        root = self.root
        nhits = self.nstates

        self.ievt +=1
        states = self.evt['gen/states']
        if (not states): return False

        pp = None
        if (self.ievt<self.nevts):
            pp = '_evt'+str(self.ievt)
        if (pp):
            root.h2d('gh_xi'+pp,nhits,0,1.*nhits,100,-15.,15.)
            root.h2d('gh_yi'+pp,nhits,0,1.*nhits,100,-15.,15.)
            root.h2d('gh_zi'+pp,nhits,0,1.*nhits,100,-5.,20.)
            root.h2d('gh_ei'+pp,nhits,0,1.*nhits,100,0.,3.)
            
            root.h2d('gh_xz'+pp,100,-5.,20,100,-15.,15.)
            root.h2d('gh_yz'+pp,100,-5.,20,100,-15.,15.)
            root.h2d('gh_ez'+pp,100,-5.,20,100,0.,3.)
            
            root.h2d('gh_kxz'+pp,100,-5.,20,100,-15.,15.)
            root.h2d('gh_kyz'+pp,100,-5.,20,100,-15.,15.)
            root.h2d('gh_kez'+pp,100,-5.,20,100,0.,3.)
        
        for i in range(len(states)-1):
            x0,y0,z0,ux0,uy0,uz0,ee0 = states[i]
            x1,y1,z1,ux1,uy1,uz1,ee1 = states[i+1]
            th0 = acos(uz0); th1 = acos(uz1);
            ph0 = 0.; ph1 = 0.;
            if (uy0 != 0.): ph0 = atan(ux0/uy0)
            if (uy1 != 0.): ph1 = atan(ux1/uy1)
            root.fill('gh_dx',x1-x0)
            root.fill('gh_dy',y1-y0)
            root.fill('gh_dz',z1-z0)
            root.fill('gh_dee',ee1-ee0)
            root.fill('gh_dth',th1-th0)
            root.fill('gh_dph',ph1-ph0)
            root.fill('gh_dxi',i+1,x1-x0)
            root.fill('gh_dyi',i+1,y1-y0)
            root.fill('gh_dzi',i+1,z1-z0)
            root.fill('gh_deei',i+1,ee1-ee0)
            root.fill('gh_dthi',i+1,th1-th0)
            root.fill('gh_dphi',i+1,ph1-ph0)

        for i,state in enumerate(states):
            x,y,z,ux,uy,uz,ee = state
            root.fill('gh_x',x)
            root.fill('gh_y',y)
            root.fill('gh_z',z)
            root.fill('gh_e',ee)
            root.fill('gh_xi',i,x)
            root.fill('gh_yi',i,y)
            root.fill('gh_zi',i,z)
            root.fill('gh_ei',i,ee)
            root.fill('gh_xz',z,x)
            root.fill('gh_yz',z,y)
            root.fill('gh_ez',z,ee)
            if (pp):
                root.fill('gh_xi'+pp,i,x)
                root.fill('gh_yi'+pp,i,y)
                root.fill('gh_zi'+pp,i,z)
                root.fill('gh_ei'+pp,i,ee)
                root.fill('gh_xz'+pp,z,x)
                root.fill('gh_yz'+pp,z,y)
                root.fill('gh_ez'+pp,z,ee)
        
        return True

class KFGenHistos(IAlg):

    def define(self):
        self.nevts = 0
        self.nhits = 30
        return
    
    def initialize(self):
        
        self.ievt = -1
        nhits = self.nhits

        self.root.h1d('gkh_nhits',nhits+1,0.,1.*(nhits+1))

        self.root.h2d('gkh_xi',nhits,0.,1.*nhits,100,-15.,15.)
        self.root.h2d('gkh_yi',nhits,0.,1.*nhits,100,-15.,15.)
        self.root.h2d('gkh_ei',nhits,0.,1.*nhits,100,0.,3.)
        
        self.root.h2d('gkh_xz',100,-5.,20,100,-15.,15.)
        self.root.h2d('gkh_yz',100,-5.,20,100,-15.,15.)
        self.root.h2d('gkh_ez',100,-5.,20,100,0.,3.)

        return True

    def execute(self):

        root = self.root
        nhits = self.nhits

        self.ievt +=1
        nodes = self.evt['gen/nodes']
        if (not nodes): return False

        pp = None
        if (self.ievt<self.nevts):
            pp = '_evt'+str(self.ievt)
        if (pp):            
            root.h2d('gkh_xi'+pp,nhits,0.,1.*nhits,100,-10,10.)
            root.h2d('gkh_yi'+pp,nhits,0.,1.*nhits,100,-10,10.)
            root.h2d('gkh_ei'+pp,nhits,0.,1.*nhits,100,0.,3.)
            root.h2d('gkh_xz'+pp,100,-5.,20,100,-15.,15.)
            root.h2d('gkh_yz'+pp,100,-5.,20,100,-15.,15.)
            root.h2d('gkh_ez'+pp,100,-5.,20,100,0.,3.)
        
        root.fill('gkh_nhits',len(nodes))

        for i,node in enumerate(nodes):
            st = node.getstate('true')
            x,y,tx,ty,ee = st.vec; z = st.zrun
            root.fill('gkh_xi',i,x)
            root.fill('gkh_yi',i,y)
            root.fill('gkh_ei',i,ee)
            root.fill('gkh_xz',z,x)
            root.fill('gkh_yz',z,y)
            root.fill('gkh_ez',z,ee)
            if (pp):
                root.fill('gkh_xi'+pp,i,x)
                root.fill('gkh_yi'+pp,i,y)
                root.fill('gkh_ei'+pp,i,ee)
                root.fill('gkh_xz'+pp,z,x)
                root.fill('gkh_yz'+pp,z,y)
                root.fill('gkh_ez'+pp,z,ee)
        return True
                    

class KFFitHistos(IAlg):

    def define(self):
        self.nhits = 20

    def initialize(self):
       
        self.root.h1d('kf_nhits',self.nhits+1,0.,self.nhits+1)       
        self.root.h1d('kf_fchi',100,0.,5.*self.nhits)
        self.root.h1d('kf_schi',100,0.,5.*self.nhits)

        self.root.h2d('kf_hxi',self.nhits,0,self.nhits,100,-10.,10.)
        self.root.h2d('kf_hyi',self.nhits,0,self.nhits,100,-10.,10.)
        self.root.h2d('kf_hzi',self.nhits,0,self.nhits,100,-5.,20.)
        self.root.h2d('kf_htxi',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('kf_htyi',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('kf_heei',self.nhits,0,self.nhits,100,0.,3.)
        self.root.h2d('kf_hfchii',self.nhits,0,self.nhits,100,0.,10.)
        self.root.h2d('kf_hschii',self.nhits,0,self.nhits,100,0.,10.)
        self.root.h2d('kf_hxpuli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('kf_hypuli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('kf_htxpuli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('kf_htypuli',self.nhits,0,self.nhits,100,-5.,5.)
        
        self.root.h1d('kf_hx',100,-10.,10.)
        self.root.h1d('kf_hy',100,-10.,10.)
        self.root.h1d('kf_hz',100,-5.,20.)
        self.root.h1d('kf_htx',100,-5.,5.)
        self.root.h1d('kf_hty',100,-5.,5.)
        self.root.h1d('kf_hee',100,0.,3.)
        self.root.h1d('kf_hfchi',100,0.,10.)
        self.root.h1d('kf_hschi',100,0.,10.)
        self.root.h1d('kf_hxpul',100,-5.,5.)
        self.root.h1d('kf_hypul',100,-5.,5.)
        self.root.h1d('kf_htxpul',100,-5.,5.)
        self.root.h1d('kf_htypul',100,-5.,5.)

    def execute(self):

        root = self.root
        kf = self.evt['kf']
        if (not kf): return False
        
        nhits = len(kf)
        fchi = kf.chi2('filter')
        schi = kf.chi2('smooth')
        self.root.fill('kf_nhits',nhits)
        self.root.fill('kf_fchi',fchi)
        self.root.fill('kf_schi',schi)

        for i in range(nhits):
            node = kf.nodes[i]
            state = node.getstate('true')
            rx,ry,rtx,rty,ree = state.vec
            rz = state.zrun
            node = kf.nodes[i]
            res = node.residual('smooth')
            x,sx = node.param('smooth',0)
            y,sy = node.param('smooth',1)
            tx,stx = node.param('smooth',2)
            ty,sty = node.param('smooth',3)
            ee,see = node.param('smooth',4)
            root.fill("kf_hxi",i,x)
            root.fill("kf_hyi",i,y)
            root.fill("kf_hzi",i,rz)
            root.fill("kf_htxi",i,tx)
            root.fill("kf_htyi",i,ty)
            root.fill("kf_heei",i,ee)
            root.fill('kf_hxpuli',i,(x-rx)/sx)
            root.fill('kf_hypuli',i,(y-ry)/sy)
            root.fill('kf_htxpuli',i,(tx-rtx)/stx)
            root.fill('kf_htypuli',i,(ty-rty)/sty)
            fchi = node.getchi2('filter')
            schi = node.getchi2('smooth')
            root.fill('kf_hfchii',i,fchi)
            root.fill('kf_hschii',i,schi)

            root.fill("kf_hx",x)
            root.fill("kf_hy",y)
            root.fill("kf_hz",rz)
            root.fill("kf_htx",tx)
            root.fill("kf_hty",ty)
            root.fill("kf_hee",ee)
            root.fill('kf_hxpul',(x-rx)/sx)
            root.fill('kf_hypul',(y-ry)/sy)
            root.fill('kf_htxpul',(tx-rtx)/stx)
            root.fill('kf_htypul',(ty-rty)/sty)
            fchi = node.getchi2('filter')
            schi = node.getchi2('smooth')
            root.fill('kf_hfchi',fchi)
            root.fill('kf_hschi',schi)
    
        return True

#-----------------------------

def ck_gen():

    alex = Alex()
    alex.nevts = 100
    root = ROOTSvc('root','ck_anextgen.root')
    alex.addsvc(root)

    radlen = 0.
    nhits = 20
    dz = 0.3
    pgas = 10.

    galg = Gen('gen')
    galg.pgas = pgas
    galg.nhits = nhits
    galg.dz = dz
    galg.xres = 1.

    ghalg = GenHistos('genhistos')  
    ghalg.imports.append('root')
    ghalg.nhits = 250

    #kgalg = KFGen('kfgen')
    #kgalg.pgas = pgas
    #kgalg.radlen = radlen
    #kgalg.nhits = nhits
    #kgalg.dz = dz

    kghalg = KFGenHistos('kfgenhistos')  
    kghalg.imports.append('root')
    kghalg.nhits = nhits

    kalg = KFFit('fit')
    kalg.pgas = pgas
    #kalg.radlen = radlen
    
    khalg = KFFitHistos('fithistos')  
    khalg.imports.append('root')
    khalg.nhits = nhits

    alex.addalg(galg)
    alex.addalg(ghalg)
    #alex.addalg(kgalg)
    alex.addalg(kghalg)
    alex.addalg(kalg)
    alex.addalg(khalg)

    alex.run()
    
    
if __name__ == '__main__':

    ck_gen()
