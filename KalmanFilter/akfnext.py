from alex.alex import IAlg
from alex.alex import Alex
from alex.rootsvc import ROOTSvc
import random
from KFBase import KFVector, KFMatrix, KFMatrixNull, KFMatrixUnitary
from kfnext import NEXT, nextgenerator, nextfilter, V0, H0
from kfnext import simplegenerator, simplefilter
from kfgenerator import zsample,zransample,zavesample,zsegments
from kffilter import randomnode, KFData
from kfzline import zstate
from math import *

"""

Alex example to check KFNextFilter

"""

class GenerateStates(IAlg):

    def define(self):
        self.pgas = 10.
        self.deltae = 0.01 # 10 keV
        self.radlen = -1
        self.E0 = 2.5
        self.ux = 0.
        self.uy = 0.
        self.uz = +1.
        self.emin = 0.5
        return

    def execute(self):
        next = NEXT(self.pgas)
        kfgen = nextgenerator(next,deltae=self.deltae,emin=self.emin)
        if (self.radlen>0):
            kfgen = simplegenerator(radlen=self.radlen,deltae=self.deltae,emin=self.emin)
        state0 = [0.,0.,0.,self.ux,self.uy,self.uz,self.E0]
        states = kfgen.generate(state0)
        self.evt['sim/states'] = states
        #for state in states:
        #    print ' state ',state
        return True

class GenerateDigits(IAlg):

    def define(self):
        self.dz = .5
        self.ZZ = 30.
        return

    def initialize(self):
        self.nhits = int(2.*self.ZZ/self.dz)
        self.zs = map(lambda i: -1.*self.ZZ+i*self.dz,range(self.nhits))
        self.depsilon = self.dz/100.
        print 'zs ',self.zs
        return True

    def execute(self):
        states = self.evt['sim/states']
        if (not states): return False
        zstates = zsample(states,self.zs)
        ok = len(zstates)>0
        if (ok): self.evt['sim/digits']=zstates
        #for izstate in zstates:
        #    print 'zstate ',izstate
        return ok


class GenerateDigits2(IAlg):

    def define(self):
        self.p0 = .2
        return

    def execute(self):
        print 'keys ',self.evt.keys()
        states = self.evt['sim/states']
        if (not states): return False
        zstates = zransample(states,self.p0)
        ok = len(zstates)>0
        if (ok): self.evt['sim/digits']=zstates
        #for izstate in zstates:
        #    print 'zstate ',izstate
        return ok

class GenerateDigits3(IAlg):

    def define(self):
        self.nave = 4
        return

    def execute(self):
        print 'keys ',self.evt.keys()
        states = self.evt['sim/states']
        if (not states): return False
        zstates = zavesample(states,self.nave)
        ok = len(zstates)>0
        if (ok): self.evt['sim/digits']=zstates
        #for izstate in zstates:
        #    print 'zstate ',izstate
        return ok


class GenerateSegments(IAlg):

    def define(self):
        self.unique = False

    def execute(self):
        states = self.evt['sim/digits']     
        if (not states): return False
        if (self.unique):
            segments = [states,]
        else:
            segments = zsegments(states)
        print " # of segments ",len(segments)
        ok = (len(segments)>=1)
        if (ok): self.evt['sim/segments'] = segments
        #for seg in segments:
        #    print " segment: "
        #    for state in seg: print state
        return ok


class GenerateNodes(IAlg):

    def define(self):
        self.xres = 0.1
        self.minnodes = 4
        return

    def initialize(self):
        self.V = V0(self.xres)
        return True

    def execute(self):
        segments = self.evt['sim/segments']   
        if (not segments): return False
        segnodes = []
        for seg in segments:
            if (len(seg)<=self.minodes): continue
            zstates = map(zstate,seg)
            nodes = map(lambda st: randomnode(st,H0,self.V),zstates)
            segnodes.append(nodes)
        ok = (len(segnodes)>0)
        if ok: self.evt['sim/segnodes'] = segnodes
        #for znode in segnodes:
        #    print ' segnode ',len(znode)
        #    print ' zruns ',map(lambda nd: nd.hit.zrun,znode)
        return ok


class GenenerateKFNodes(IAlg):
    """ Algorithm to generate SL tracks, fit them with KFNextFilter and fill histos
    """

    def define(self):
        self.pgas = 10.
        self.dz = 0.5
        self.ZZ = 30.
        self.xres = 0.1
        self.E0 = 2.5
        return

    def initialize(self):
        self.nhits = int(2.*self.ZZ/self.dz)
        self.zs = map(lambda i: -1.*self.ZZ+i*self.dz,range(self.nhits))
        self.V = V0(self.xres)  
        return True

    def execute(self):
        # Generate Track
        #----------------
        #empty this
        nullhits = map(lambda z: KFData(KFVector([0,0]),self.V,z),self.zs)
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
            self.evt['sim/nodes']=knodes
        return ok


class KFFit(IAlg):
    """ Algorithm to generate SL tracks, fit them with KFNextFilter and fill histos
    """

    def define(self):
        self.pgas = 10.
        self.radlen = -1
        self.E0 = 2.5
        return


    def fitsegment(self,nodes):
        if (len(nodes)<=2): return 

        next = NEXT(self.pgas)
        kf = nextfilter(next)
        if (self.radlen>0.):
            self.msg.info('KF Simple radlen',self.radlen)
            kf = simplefilter(self.radlen)
        #kf = simplefilter(self.radlen)
        kf.clear()
        kf.setnodes(nodes)
        print " fitsegment nodes ",len(kf.nodes) 

        # fit!
        x0true = nodes[0].getstate('true').vec
        x0 = KFVector(x0true)
        z0 = nodes[0].zrun
        dz = nodes[1].zrun-z0
        if (dz==0.): 
            print ' Fit failed, not valid input nodes ',z0,dz
            return False,None
        zdir = dz/abs(dz)
        C1 = 1.*KFMatrixUnitary(5)
        state0 = KFData(x0,C1,z0-0.1*zdir,pars={'uz':zdir})
        ok,fchi,schi = kf.fit(state0)  
        if (not ok):
            print " Fit Failed! "
            return False,None
        print ' Fit chi2 ',dz,fchi,schi
        return ok,kf

    def execute(self):
        segnodes = self.evt['sim/segnodes']
        if (not segnodes): return False
        kfs = []
        for nodes in segnodes:
            print ' nodes to fit ',len(nodes)
            ok,kf = self.fitsegment(nodes)
            if (ok): kfs.append(kf)
        ok = (len(kfs)>0)
        if (ok): self.evt['rec/kfs']=kfs
        print " segments fits ",len(kfs)
        return ok

#--------- histogramming!!

class HistosGenerateStates(IAlg):

    def define(self):
        self.nevts = 0
        self.nstates = 250
        self.prefix = 'gh_'
        return
    
    def initialize(self):
        
        self.ievt = -1
        nhits = self.nstates

        self.root.h1d(self.prefix+'x',100,-15.,15.)
        self.root.h1d(self.prefix+'y',100,-15.,15.)
        self.root.h1d(self.prefix+'z',100,-10.,20.)
        self.root.h1d(self.prefix+'e',100,0.,3.)

        self.root.h2d(self.prefix+'xi',nhits,0,1.*nhits,100,-15.,15.)
        self.root.h2d(self.prefix+'yi',nhits,0,1.*nhits,100,-15.,15.)
        self.root.h2d(self.prefix+'zi',nhits,0,1.*nhits,100,-5.,20.)
        self.root.h2d(self.prefix+'ei',nhits,0,1.*nhits,100,0.,3.)

        self.root.h1d(self.prefix+'dx',100,-0.5,0.5)
        self.root.h1d(self.prefix+'dy',100,-0.5,0.5)
        self.root.h1d(self.prefix+'dz',100,-0.5,0.5)
        self.root.h1d(self.prefix+'dth',100,-pi,pi)
        self.root.h1d(self.prefix+'dph',100,-pi,pi)
        self.root.h1d(self.prefix+'dee',100,-0.2,0.2)

        self.root.h2d(self.prefix+'dxi',nhits,0.,1.*nhits,100,-0.5,0.5)
        self.root.h2d(self.prefix+'dyi',nhits,0.,1.*nhits,100,-0.5,0.5)
        self.root.h2d(self.prefix+'dzi',nhits,0.,1.*nhits,100,-0.5,0.5)
        self.root.h2d(self.prefix+'dthi',nhits,0.,1.*nhits,100,-pi,pi)
        self.root.h2d(self.prefix+'dphi',nhits,0.,1.*nhits,100,-pi,pi)
        self.root.h2d(self.prefix+'deei',nhits,0.,1.*nhits,100,-0.2,0.2)
        
        self.root.h2d(self.prefix+'xz',100,-5.,20,100,-15.,15.)
        self.root.h2d(self.prefix+'yz',100,-5.,20,100,-15.,15.)
        self.root.h2d(self.prefix+'ez',100,-5.,20,100,0.,3.)

    def execute(self):

        root = self.root
        nhits = self.nstates

        self.ievt +=1
        states = self.evt['sim/states']
        if (not states): return False

        pp = None
        if (self.ievt<self.nevts):
            pp = '_evt'+str(self.ievt)
        if (pp):
            root.h2d(self.prefix+'xi'+pp,nhits,0,1.*nhits,100,-15.,15.)
            root.h2d(self.prefix+'yi'+pp,nhits,0,1.*nhits,100,-15.,15.)
            root.h2d(self.prefix+'zi'+pp,nhits,0,1.*nhits,100,-5.,20.)
            root.h2d(self.prefix+'ei'+pp,nhits,0,1.*nhits,100,0.,3.)
            
            root.h2d(self.prefix+'xz'+pp,100,-5.,20,100,-15.,15.)
            root.h2d(self.prefix+'yz'+pp,100,-5.,20,100,-15.,15.)
            root.h2d(self.prefix+'ez'+pp,100,-5.,20,100,0.,3.)
            
            root.h2d(self.prefix+'kxz'+pp,100,-5.,20,100,-15.,15.)
            root.h2d(self.prefix+'kyz'+pp,100,-5.,20,100,-15.,15.)
            root.h2d(self.prefix+'kez'+pp,100,-5.,20,100,0.,3.)
        
        for i in range(len(states)-1):
            x0,y0,z0,ux0,uy0,uz0,ee0 = states[i]
            x1,y1,z1,ux1,uy1,uz1,ee1 = states[i+1]
            th0 = acos(uz0); th1 = acos(uz1);
            ph0 = 0.; ph1 = 0.;
            if (uy0 != 0.): ph0 = atan(ux0/uy0)
            if (uy1 != 0.): ph1 = atan(ux1/uy1)
            root.fill(self.prefix+'dx',x1-x0)
            root.fill(self.prefix+'dy',y1-y0)
            root.fill(self.prefix+'dz',z1-z0)
            root.fill(self.prefix+'dee',ee1-ee0)
            root.fill(self.prefix+'dth',th1-th0)
            root.fill(self.prefix+'dph',ph1-ph0)
            root.fill(self.prefix+'dxi',i+1,x1-x0)
            root.fill(self.prefix+'dyi',i+1,y1-y0)
            root.fill(self.prefix+'dzi',i+1,z1-z0)
            root.fill(self.prefix+'deei',i+1,ee1-ee0)
            root.fill(self.prefix+'dthi',i+1,th1-th0)
            root.fill(self.prefix+'dphi',i+1,ph1-ph0)

        for i,state in enumerate(states):
            x,y,z,ux,uy,uz,ee = state
            root.fill(self.prefix+'x',x)
            root.fill(self.prefix+'y',y)
            root.fill(self.prefix+'z',z)
            root.fill(self.prefix+'e',ee)
            root.fill(self.prefix+'xi',i,x)
            root.fill(self.prefix+'yi',i,y)
            root.fill(self.prefix+'zi',i,z)
            root.fill(self.prefix+'ei',i,ee)
            root.fill(self.prefix+'xz',z,x)
            root.fill(self.prefix+'yz',z,y)
            root.fill(self.prefix+'ez',z,ee)
            if (pp):
                root.fill(self.prefix+'xi'+pp,i,x)
                root.fill(self.prefix+'yi'+pp,i,y)
                root.fill(self.prefix+'zi'+pp,i,z)
                root.fill(self.prefix+'ei'+pp,i,ee)
                root.fill(self.prefix+'xz'+pp,z,x)
                root.fill(self.prefix+'yz'+pp,z,y)
                root.fill(self.prefix+'ez'+pp,z,ee)
        
        return True

class HistosGenerateSegments(IAlg):

    def initialize(self):
        self.root.h1d('gs_nseg',20,0.,20.)
        self.root.h1d('gs_segsize',100,0.,100.)
        self.root.h1d('gs_ei',100,0.,3.)
        self.root.h2d('gs_sizeei',100,0.,3.,100,0.,100.)
        return True  

    def execute(self):
        segments = self.evt['sim/segments']
        if (not segments): return False
        self.root.fill('gs_nseg',len(segments))
        for seg in segments:
            self.root.fill('gs_segsize',len(seg))
            self.root.fill('gs_ei',seg[0][-1])
            self.root.fill('gs_sizeei',seg[0][-1],len(seg))
        return True


class HistosGenerateNodes(IAlg):

    def define(self):
        self.nevts = 0
        self.nhits = 30
        self.prefix = 'gn_'
        return
    
    def initialize(self):
        
        self.ievt = -1
        nhits = self.nhits

        self.root.h1d(self.prefix+'nhits',nhits+1,0.,1.*(nhits+1))

        self.root.h2d(self.prefix+'xi',nhits,0.,1.*nhits,100,-15.,15.)
        self.root.h2d(self.prefix+'yi',nhits,0.,1.*nhits,100,-15.,15.)
        self.root.h2d(self.prefix+'ei',nhits,0.,1.*nhits,100,0.,3.)
        
        self.root.h2d(self.prefix+'xz',100,-5.,20,100,-15.,15.)
        self.root.h2d(self.prefix+'yz',100,-5.,20,100,-15.,15.)
        self.root.h2d(self.prefix+'ez',100,-5.,20,100,0.,3.)

        return True

    def execute(self):

        root = self.root
        nhits = self.nhits

        self.ievt +=1
        nodes = self.evt['sim/nodes']
        if (not nodes): return False

        pp = None
        if (self.ievt<self.nevts):
            pp = '_evt'+str(self.ievt)
        if (pp):            
            root.h2d(self.prefix+'xi'+pp,nhits,0.,1.*nhits,100,-10,10.)
            root.h2d(self.prefix+'yi'+pp,nhits,0.,1.*nhits,100,-10,10.)
            root.h2d(self.prefix+'ei'+pp,nhits,0.,1.*nhits,100,0.,3.)
            root.h2d(self.prefix+'xz'+pp,100,-5.,20,100,-15.,15.)
            root.h2d(self.prefix+'yz'+pp,100,-5.,20,100,-15.,15.)
            root.h2d(self.prefix+'ez'+pp,100,-5.,20,100,0.,3.)
        
        root.fill(self.prefix+'nhits',len(nodes))

        for i,node in enumerate(nodes):
            st = node.getstate('true')
            x,y,tx,ty,ee = st.vec; z = st.zrun
            root.fill(self.prefix+'xi',i,x)
            root.fill(self.prefix+'yi',i,y)
            root.fill(self.prefix+'ei',i,ee)
            root.fill(self.prefix+'xz',z,x)
            root.fill(self.prefix+'yz',z,y)
            root.fill(self.prefix+'ez',z,ee)
            if (pp):
                root.fill(self.prefix+'xi'+pp,i,x)
                root.fill(self.prefix+'yi'+pp,i,y)
                root.fill(self.prefix+'ei'+pp,i,ee)
                root.fill(self.prefix+'xz'+pp,z,x)
                root.fill(self.prefix+'yz'+pp,z,y)
                root.fill(self.prefix+'ez'+pp,z,ee)
        return True

                  
class HistosKFFit(IAlg):

    def define(self):
        self.nhits = 40
        self.prefix = 'kf_'
        self.emin = 0.5
        self.forward = True
        self.backward = True
        return

    def initialize(self):
       
        self.root.h1d(self.prefix+'nfits',10,0.,10)

        self.root.h1d(self.prefix+'nhits',self.nhits+1,0.,self.nhits+1)       
        self.root.h1d(self.prefix+'fchi',100,0.,5.*self.nhits)
        self.root.h1d(self.prefix+'schi',100,0.,5.*self.nhits)

        self.root.h2d(self.prefix+'hxi',self.nhits,0,self.nhits,100,-10.,10.)
        self.root.h2d(self.prefix+'hyi',self.nhits,0,self.nhits,100,-10.,10.)
        self.root.h2d(self.prefix+'hzi',self.nhits,0,self.nhits,100,-5.,20.)
        self.root.h2d(self.prefix+'htxi',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'htyi',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'heei',self.nhits,0,self.nhits,100,0.,3.)
        self.root.h2d(self.prefix+'hfchii',self.nhits,0,self.nhits,100,0.,10.)
        self.root.h2d(self.prefix+'hschii',self.nhits,0,self.nhits,100,0.,10.)
        self.root.h2d(self.prefix+'hxpuli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'hypuli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'htxpuli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'htypuli',self.nhits,0,self.nhits,100,-5.,5.)
        
        self.root.h1d(self.prefix+'hx',100,-10.,10.)
        self.root.h1d(self.prefix+'hy',100,-10.,10.)
        self.root.h1d(self.prefix+'hz',100,-5.,20.)
        self.root.h1d(self.prefix+'htx',100,-5.,5.)
        self.root.h1d(self.prefix+'hty',100,-5.,5.)
        self.root.h1d(self.prefix+'hee',100,0.,3.)
        self.root.h1d(self.prefix+'hfchi',100,0.,10.)
        self.root.h1d(self.prefix+'hschi',100,0.,10.)
        self.root.h1d(self.prefix+'hxpul',100,-5.,5.)
        self.root.h1d(self.prefix+'hypul',100,-5.,5.)
        self.root.h1d(self.prefix+'htxpul',100,-5.,5.)
        self.root.h1d(self.prefix+'htypul',100,-5.,5.)

    def execute(self):
        kfs = self.evt['rec/kfs']
        if (not kfs): return False
        self.root.fill(self.prefix+'nfits',len(kfs))
        for kf in kfs: self.kffill(kf)
        return True

    def kffill(self,kf):
        root = self.root
        nhits = len(kf)
        dz = kf.nodes[0].getstate('true').pars['uz']
        if (dz>0 and not self.forward): return
        if (dz<0 and not self.backward): return
        ene = kf.nodes[0].getstate('true').vec[4]
        if (ene<self.emin): return
        fchi = kf.chi2('filter')
        schi = kf.chi2('smooth')
        self.root.fill(self.prefix+'nhits',nhits)
        self.root.fill(self.prefix+'fchi',fchi)
        self.root.fill(self.prefix+'schi',schi)

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
            root.fill(self.prefix+'hxpuli',i,(x-rx)/sx)
            root.fill(self.prefix+'hypuli',i,(y-ry)/sy)
            root.fill(self.prefix+'htxpuli',i,(tx-rtx)/stx)
            root.fill(self.prefix+'htypuli',i,(ty-rty)/sty)
            fchi = node.getchi2('filter')
            schi = node.getchi2('smooth')
            root.fill(self.prefix+'hfchii',i,fchi)
            root.fill(self.prefix+'hschii',i,schi)

            root.fill("kf_hx",x)
            root.fill("kf_hy",y)
            root.fill("kf_hz",rz)
            root.fill("kf_htx",tx)
            root.fill("kf_hty",ty)
            root.fill("kf_hee",ee)
            root.fill(self.prefix+'hxpul',(x-rx)/sx)
            root.fill(self.prefix+'hypul',(y-ry)/sy)
            root.fill(self.prefix+'htxpul',(tx-rtx)/stx)
            root.fill(self.prefix+'htypul',(ty-rty)/sty)
            fchi = node.getchi2('filter')
            schi = node.getchi2('smooth')
            root.fill(self.prefix+'hfchi',fchi)
            root.fill(self.prefix+'hschi',schi)
    
        return True

#-----------------------------
# Alex
#-----------------------------

def ck_gen():

    alex = Alex()
    alex.nevts = 100
    root = ROOTSvc('root','akfnext.root')
    alex.addsvc(root)

    #radlen = 1500 # if radlen>0 use simple generator and filter
    emin = 0.5 # minimum energy to stop generation

    agenstates = GenerateStates('agenstates')
    agenstates.uz = 1.
    agenstates.emin = emin
    agenstates.deltae = 0.01
    #agenstates.radlen = radlen
    alex.addalg(agenstates)

    hisgenstates = HistosGenerateStates('hisgenstates')  
    hisgenstates.imports.append('root')
    alex.addalg(hisgenstates)

    agendigits = GenerateDigits('agendigits')
    agendigits.dz = 0.5
    #agendigits.p0 = 0.25
    #agendigits.nave = 10
    alex.addalg(agendigits)  

    #hisgendigits = HistosGenerateDigits('hisgendigits')     
    #hisgen.prefix = 'hd_'
    #hisgendigits.imports.append('root')
    #alex.addalg(hisgendigits)

    agensegments = GenerateSegments('agensegments')
    agensegments.unique = False
    alex.addalg(agensegments)

    hisgensegments = HistosGenerateSegments('hisgensegments')
    hisgensegments.imports.append('root')
    alex.addalg(hisgensegments)

    agennodes = GenerateNodes('agennodes')
    agennodes.minodes = 10
    alex.addalg(agennodes)

    #hisgennodes = HistosGenerateNodes('hisgennodes')
    #hisgennodes.imports.append('root')
    #alex.addalg(hisgennodes)

    akf = KFFit('akf')
    #akf.radlen = radlen
    alex.addalg(akf)

    hisakf = HistosKFFit('hisakf')
    hisakf.imports.append('root')
    hisakf.forward = True
    hisakf.backward = True
    hisakf.emin = .5
    alex.addalg(hisakf)

    alex.run()    
    
if __name__ == '__main__':

    ck_gen()
