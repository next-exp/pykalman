from math import *
import random
from copy import deepcopy

from alex.alex import IAlg
from alex.alex import Alex
from alex.rootsvc import ROOTSvc

from KFBase import KFVector, KFMatrix, KFMatrixNull, KFMatrixUnitary
from kfnext import NEXT, nextgenerator, nextfilter, V0, H0
from kfgenerator import zavesample,zrunsample
from kffilter import randomnode, KFData, KFNode
from kfzline import zstate
#from troot import tcanvas

"""
Alex example to check KFNextFilter
"""
#------- utilities

gpgas = 10.
gnext = NEXT(gpgas)
gnextgen = nextgenerator(gnext)
gnextfit = nextfilter()

def genele(ene = 2.5):
    theta = random.uniform(0.,pi)
    phi = random.uniform(0,2.*pi)
    ux,uy,uz = sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)
    #ux,uy,uz=0.,0.,1.
    state0 = (0.,0.,0.,ux,uy,uz,ene)
    states = gnextgen.generate(state0)
    # print 'genele-states-'states
    return states

def genbeta(ene = 2.5):
    emin = 0.512
    ei0 = random.uniform(emin,ene-emin)
    ei1 = ene-ei0
    e0 = min(ei0,ei1)
    e1 = max(ei0,ei1)
    states1 = genele(e0)
    states2 = genele(e1)
    # print 'genbeta-states1 ',states1
    # print 'genbeta-states2 ',states2
    return [states1,states2]

def createzdigits(states,zs):
    #zstates = zrunsample(states,zs)
    zstates = zavesample(states)
    denes = map(lambda i: zstates[i][-1]-zstates[i+1][-1],range(len(zstates)-1))    
    denes.append(zstates[-1][-1])
    digits = map(lambda st,dene: [st[0],st[1],st[2],dene],zstates,denes)
    # print 'createzdigits digits'
    return (digits,zstates)

def createhits(digits,xres=0.1):
    hits = []
    for digit in digits:
        x,y,z,dene = digit
        x = x+random.gauss(0.,xres)
        y = y+random.gauss(0.,xres)
        V = (xres*xres)*KFMatrixUnitary(2)
        hit = KFData(KFVector([x,y]),V,z)
        hit.dene = dene
        hits.append(hit)
    # print 'createhits ',hits
    return hits

def createnodes(hits,H=H0,states=None):
    nodes = map(lambda hit: KFNode(hit,H),hits)
    if (states):
        zstates = map(zstate,states)
        for i,node in enumerate(nodes):
            node.setstate('true',zstates[i])
    #print "create nodes "
    return nodes

def preparenodes(nodes):
    denes = map(lambda nd: nd.hit.dene,nodes)
    for i,node in enumerate(nodes):
        node.hit.ene = sum(denes[i:])
    enes = map(lambda nd: nd.hit.ene,nodes)
    #print ' preparenodes enes ',enes
    return

def seedstate(nodes):
    x0 = nodes[0].hit.vec
    x1 = nodes[1].hit.vec
    z0 = nodes[0].zrun
    dz = nodes[1].zrun-z0
    if (dz==0.): 
        print ' Fit failed, not valid input nodes ',z0,dz
        return False,None
    ene = nodes[0].hit.ene
    tx,ty = (x0[0]-x1[0])/dz,(x0[1]-x1[1])/dz
    xx = KFVector([x0[0],x0[1],tx,ty,ene])
    zdir = dz/abs(dz)
    C1 = 10.*KFMatrixUnitary(5)
    state0 = KFData(xx,C1,z0-0.1*zdir,pars={'uz':zdir})
    #print "seedstate ",state0
    return state0

def fitnodes(nodes):
    preparenodes(nodes)
    state0 = seedstate(nodes)
    gnextfit.clear()
    gnextfit.setnodes(nodes)
    cc = gnextfit.fit(state0)
    kf = deepcopy(gnextfit)
    #print 'fitnodes ',cc    
    return cc,kf

def chiblocks(nodes,maxchi=5,mtype='filter',nblock=5):
    chis = map(lambda nd: nd.getchi2(mtype),nodes)
    chis = filter(lambda chi: chi<maxchi,chis)
    nn = len(chis)
    nb = int(nn/(1.*nblock))
    bchis = map(lambda i: sum(chis[i*nblock:(i+1)*nblock])/(2.*nblock-4.),range(nb))
    print 'chiblocks ',nn,nb,bchis
    #chis.append(sum(chi2[nblock*nb:])/(2*nblock-4))
    return bchis

def chindf(nodes,maxchi=30.,mtype='filter'):
    chi = map(lambda nd: nd.getchi2(mtype),nodes)
    chi = filter(lambda ch: ch<maxchi,chi)
    nn = len(chi)
    if nn<=2: chi = 1.
    else: chi = sum(chi)/(2.*nn-4.)
    return chi

def achi(nodes,maxchi=30.,m=3):
    nn = len(nodes)
    ni = int(nn/m)
    if (ni<=2): return 1.,1.
    acs = []
    for mtype in ['filter','smooth']:
        ch1 = chindf(nodes[:ni],mtype=mtype)
        ch3 = chindf(nodes[(m-1)*ni:],mtype=mtype)
        ac = (ch1-ch3)/(ch1+ch3)
        acs.append(ac)
    # print "achi ",acs
    return acs

def acvertex(nodes,rnodes):
    if (len(nodes)!=len(rnodes)):
        print " WARNING! acvertex, different # nodes ",len(nodes),len(rnodes)
    nn = len(nodes)
    i0,ie = 0,nn
    facs,racs = [],[]
    for i in range(i0,ie):
        fac = achi(nodes[i:ie])[0]
        rac = achi(rnodes[nn-i:ie])[0]
        facs.append(fac)
        racs.append(rac)
    print ' facs ',facs
    print ' racs ',racs
    return facs,racs

#------- algorithms
class GenerateBeta(IAlg):
    """ Algorithm to generate states (x,y,z,ux,uy,uz,p)
    """
    def define(self):
        self.E0 = 2.5 # MeV
        return

    def execute(self):
        states = genele(self.E0)
        self.evt['sim/states'] = [states,]
        val = map(lambda st: (len(st),st[0],st[-1]),states)
        self.msg.verbose(self.name,val)
        return True

class GenerateDoubleBeta(IAlg):

    def define(self):
        self.E0 = 2.5 # MeV
        return

    def execute(self):
        states = genbeta(self.E0)
        self.evt['sim/states'] = states
        data = map(lambda seg: (len(seg),seg[0],seg[-1]),states)
        self.msg.info(self.name,' states ',data)
        return True

class CreateDigits(IAlg):

    def define(self):
        self.dz = 0.5
        self.ZZ = 30.
        return

    def initialize(self):
        n = int(2.*self.ZZ/self.dz)
        self.zs = map(lambda i: -1.*self.ZZ+i*self.dz,range(n))
        self.msg.info(self.name,' zs ',self.zs)
        return True

    def execute(self):
        states = self.evt['sim/states']
        if (not states): return False
        digs,zsts = [],[]
        for seg in states:
            dig,zst = createzdigits(seg,self.zs)
            digs.append(dig)
            zsts.append(zst)
        ok = len(digs)>0
        if (ok): 
            self.evt['sim/digits']=digs
            self.evt['sim/zstates']=zsts
        data1 = map(lambda dig: (len(dig),dig[0],dig[-1]),digs)
        data2 = map(lambda zst: (len(zst),zst[0],zst[-1]),zsts)
        self.msg.info(self.name,' digits ',data1)
        self.msg.info(self.name,' zstats ',data2)    
        return True


class CreateNodes(IAlg):

    def define(self):
        self.xres = 0.1
        return

    def execute(self):
        digs = self.evt['sim/digits']
        zsts = self.evt['sim/zstates']
        if (not digs or not zsts): return False
        hits, nods = [],[]
        for i,dig in enumerate(digs):
            zst = zsts[i]
            hit = createhits(dig)
            nod = createnodes(hit,H0,zst)
            hits.append(hit)
            nods.append(nod)
        ok = len(nods)>0
        if (ok):
            self.evt['rec/hits'] = hits
            self.evt['rec/nodes'] = nods
        dat1 = map(lambda hit: (len(hit),hit[0]),hits)
        #dat2 = map(lambda nod: (len(nod),nod[0]),nods)
        self.msg.info(self.name,'hits ',dat1)
        #self.msg.info(self.name,'nodes ',dat2)
        return ok

class ReverseNodes(IAlg):

    def define(self):
        self.path = 'rec/nodes'
        self.opath = 'rec/nodes/rev'
        return True

    def execute(self):
        nodes = self.evt[self.path]
        if (not nodes): return False
        onodes = []
        for node in nodes:
            onode = deepcopy(node)
            onode.reverse()
            onodes.append(onode)
        self.evt[self.opath] = onodes
        val = map(lambda nd: (len(nd),nd[0].zrun),onodes)
        self.msg.info(self.name,' nodes ',val)
        return True

class DoubleBetaNodes(IAlg):

    def execute(self):
        digs = self.evt['sim/digits']
        zsts = self.evt['sim/zstates']
        if (not digs or not zsts): return False
        print 'digs ',len(digs),len(digs[0]),len(digs[1])
        if (len(digs)!=2): return False
        dig = deepcopy(digs[0])
        zst = deepcopy(zsts[0])
        dig.reverse(); zst.reverse()
        dig+=deepcopy(digs[1]); zst+=deepcopy(zsts[1])
        print ' digs - total ',len(dig)
        hits = [createhits(dig),]
        nods = [createnodes(hits[0],H0,zst),]
        ok = len(nods)>0
        if (ok):
            self.evt['rec/hits/bb'] = hits
            self.evt['rec/nodes/bb'] = nods
        dat1 = map(lambda hit: (len(hit),hit[0]),hits)
        #dat2 = map(lambda nod: (len(nod),nod[0]),nods)
        self.msg.info(self.name,'hits ',dat1)
        #self.msg.info(self.name,'nodes ',dat2)
        return ok

class FitNodes(IAlg):

    def define(self):
        self.path = 'rec/nodes'
        self.opath = 'rec/kfs'
        return True

    def execute(self):
        segs = self.evt[self.path]
        if (not segs): return False
        kfs = []
        for seg in segs:
            cc,kf = fitnodes(seg)
            ok,fchi,schi = cc
            self.msg.info(self.name,' fit ',cc)
            if ok: kfs.append(kf)
        ok = len(kfs)>0
        if (ok): self.evt[self.opath]=kfs
        vals = map(lambda kf: (kf.cleanchi2('filter'),kf.chi2('smooth')),kfs)
        self.msg.info(self.name,' chi2 ',vals)
        return ok

class FitScanNodes(IAlg):

    def define(self):
        self.path = 'rec/nodes/bb'
        self.opath = 'rec/kfs/scan'
        self.iscan = 2
        return

    def execute(self):
        segs = self.evt[self.path]
        digs = self.evt['sim/digits']
        if (not segs or not digs): return False
        nver = len(digs[0])
        nodes = segs[0]
        ntot = len(nodes)
        inodes = map(lambda i: self.iscan*i,range(1,int(ntot/self.iscan)-1))
        if (not (nver in inodes)): inodes.append(nver)
        inodes.sort()
        inodes = filter(lambda i: abs(i-nver)<10,inodes)
        print 'nnodes nvertex',ntot,nver
        print 'inodes ',inodes
        vals = {}
        for inode in inodes:
            rnodes = deepcopy(nodes[:inode])
            rnodes.reverse()
            fnodes = deepcopy(nodes[inode:])
            fcc,fkf = fitnodes(fnodes)
            rcc,rkf = fitnodes(rnodes)
            vals[inode-nver]=(rkf,fkf)
            self.msg.info('scan i ',inode-nver)
            self.msg.info('scan r ',rcc,achi(rkf.nodes))
            self.msg.info('scan f ',fcc,achi(fkf.nodes))
        self.evt['rec/kfs/scan']=vals
        return True

#--------- histogramming!!

class HistosGenerateStates(IAlg):

    def define(self):
        self.nstates = 250 # max number of states
        self.prefix = 'gsta_'
        return
    
    def initialize(self):
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

    def histoseg(self,states):
        root = self.root
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
        return True
    
    def execute(self):

        segs = self.evt['sim/states']
        if (not segs): return False
        for seg in segs:
            self.histoseg(seg)
        return True

                  
class HistosKFFit(IAlg):

    def define(self):
        self.nhits = 40 # max number of hits on track
        self.maxchi2 = 30.
        self.path = 'rec/kfs'
        self.prefix = 'kf_'
        self.ftype = 'filter'
        return

    def initialize(self):
       
        self.root.h1d(self.prefix+'nfits',10,0.,10)

        self.root.h1d(self.prefix+'nhits',self.nhits+1,0.,self.nhits+1)       
        self.root.h1d(self.prefix+'fchi2',100,0.,5.*self.nhits)
        self.root.h1d(self.prefix+'schi2',100,0.,5.*self.nhits)
        self.root.h1d(self.prefix+'nfhits',self.nhits,0.,self.nhits)
        self.root.h1d(self.prefix+'nshits',self.nhits,0.,self.nhits)
        self.root.h1d(self.prefix+'nfrejhits',self.nhits,0.,self.nhits)
        self.root.h1d(self.prefix+'nsrejhits',self.nhits,0.,self.nhits)
        self.root.h1d(self.prefix+'fchi2ndf',100,0.,10.)
        self.root.h1d(self.prefix+'schi2ndf',100,0.,10.)
        self.root.h1d(self.prefix+'fac',100,-1.,1.)
        self.root.h1d(self.prefix+'sac',100,-1.,1.)

        self.root.h1d(self.prefix+'x',100,-30.,30.)
        self.root.h1d(self.prefix+'y',100,-30.,30.)
        self.root.h1d(self.prefix+'z',100,-30.,30.)
        self.root.h1d(self.prefix+'tx',100,-5.,5.)
        self.root.h1d(self.prefix+'ty',100,-5.,5.)
        self.root.h1d(self.prefix+'ee',100,0.,3.)
        self.root.h1d(self.prefix+'fchi',100,0.,10.)
        self.root.h1d(self.prefix+'schi',100,0.,10.)
        self.root.h1d(self.prefix+'xpull',100,-5.,5.)
        self.root.h1d(self.prefix+'ypull',100,-5.,5.)
        self.root.h1d(self.prefix+'txpull',100,-5.,5.)
        self.root.h1d(self.prefix+'typull',100,-5.,5.)
        self.root.h1d(self.prefix+'nhits',100,-5.,5.)

        self.root.h2d(self.prefix+'xi',self.nhits,0,self.nhits,100,-10.,10.)
        self.root.h2d(self.prefix+'yi',self.nhits,0,self.nhits,100,-10.,10.)
        self.root.h2d(self.prefix+'zi',self.nhits,0,self.nhits,100,-5.,20.)
        self.root.h2d(self.prefix+'txi',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'tyi',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'eei',self.nhits,0,self.nhits,100,0.,3.)
        self.root.h2d(self.prefix+'fchii',self.nhits,0,self.nhits,100,0.,10.)
        self.root.h2d(self.prefix+'schii',self.nhits,0,self.nhits,100,0.,10.)
        self.root.hprf(self.prefix+'fchii_pf',self.nhits,0,self.nhits,0.,40.)
        self.root.hprf(self.prefix+'schii_pf',self.nhits,0,self.nhits,0.,40.)
        self.root.h2d(self.prefix+'xpulli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'ypulli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'txpulli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'typulli',self.nhits,0,self.nhits,100,-5.,5.)

    def execute(self):
        kfs = self.evt[self.path]
        if (not kfs): return False
        self.root.fill(self.prefix+'nfits',len(kfs))
        for kf in kfs: self.kffill(kf)
        return True

    def kffill(self,kf):
        root = self.root
        nhits = len(kf)
        nf,fchi = kf.cleanchi2('filter',self.maxchi2)
        ns,schi = kf.cleanchi2('smooth',self.maxchi2)
        self.root.fill(self.prefix+'nhits',nhits)
        self.root.fill(self.prefix+'fchi2',fchi)
        self.root.fill(self.prefix+'schi2',schi)
        self.root.fill(self.prefix+'nfhits',nf)
        self.root.fill(self.prefix+'nsrejhits',nhits-ns)
        self.root.fill(self.prefix+'nfrejhits',nhits-nf)
        self.root.fill(self.prefix+'nshits',ns)
        self.root.fill(self.prefix+'schi2ndf',schi/(2.*ns-4.))
        self.root.fill(self.prefix+'fchi2ndf',fchi/(2.*nf-4.))
        fac,sac = achi(kf.nodes)
        self.root.fill(self.prefix+'fac',fac)
        self.root.fill(self.prefix+'sac',sac)

        for i in range(nhits):
            node = kf.nodes[i]
            state = node.getstate('true')
            rx,ry,rtx,rty,ree = state.vec
            rz = state.zrun
            node = kf.nodes[i]
            fchi = node.getchi2('filter')
            schi = node.getchi2('smooth')
            res = node.residual(self.ftype)
            x,sx = node.param(self.ftype,0)
            y,sy = node.param(self.ftype,1)
            tx,stx = node.param(self.ftype,2)
            ty,sty = node.param(self.ftype,3)
            ee,see = node.param(self.ftype,4)

            root.fill(self.prefix+"x",x)
            root.fill(self.prefix+"y",y)
            root.fill(self.prefix+'z',rz)
            root.fill(self.prefix+"tx",tx)
            root.fill(self.prefix+"ty",ty)
            root.fill(self.prefix+"ee",ee)
            if (fchi<=self.maxchi2):
                root.fill(self.prefix+'fchi',fchi)
            if (schi<=self.maxchi2):
                root.fill(self.prefix+'schi',schi)
                root.fill(self.prefix+'xpull',(x-rx)/sx)
                root.fill(self.prefix+'ypull',(y-ry)/sy)
                root.fill(self.prefix+'txpull',(tx-rtx)/stx)
                root.fill(self.prefix+'typull',(ty-rty)/sty)

            root.fill(self.prefix+"xi",i,x)
            root.fill(self.prefix+"yi",i,y)
            root.fill(self.prefix+"zi",i,rz)
            root.fill(self.prefix+"txi",i,tx)
            root.fill(self.prefix+"tyi",i,ty)
            root.fill(self.prefix+"eei",i,ee)
            if (schi<=self.maxchi2):
                root.fill(self.prefix+'xpulli',i,(x-rx)/sx)
                root.fill(self.prefix+'ypulli',i,(y-ry)/sy)
                root.fill(self.prefix+'txpulli',i,(tx-rtx)/stx)
                root.fill(self.prefix+'typulli',i,(ty-rty)/sty)
                root.fill(self.prefix+'schii',i,schi)
                root.fill(self.prefix+'schii_pf',i,schi)
            if (fchi<=self.maxchi2):
                root.fill(self.prefix+'fchii',i,fchi)
                root.fill(self.prefix+'fchii_pf',i,fchi)
        
        return True

class AnaReverse(IAlg):

    def define(self):
        self.prefix = 'nr_'
        self.pathnor = 'rec/kfs'
        self.pathrev = 'rec/kfs/rev'
        return

    def initialize(self):
        self.root.h2d(self.prefix+'fac',40,-1.,1.,40,-1.,1.)
        self.root.h2d(self.prefix+'sac',40,-1.,1.,40,-1.,1.)
        return True

    def execute(self):
        nkfs = self.evt[self.pathnor]
        rkfs = self.evt[self.pathrev]
        if (not nkfs or not rkfs): return False
        pars = zip(nkfs,rkfs)
        for par in pars:
            nkf,rkf = par
            nac = achi(nkf.nodes)
            rac = achi(rkf.nodes)
            self.root.fill(self.prefix+'fac',nac[0],rac[0])
            self.root.fill(self.prefix+'sac',nac[1],rac[1])
        return True

class HistosScanNodes(IAlg):

    def define(self):
        self.prefix = 'vtx_'
        self.path = 'rec/kfs/scan'
        self.inod = 10
        return

    def initialize(self):

        self.root.h2d(self.prefix+'faci',60,-self.inod,self.inod,40,-2.,2.)
        self.root.h2d(self.prefix+'raci',60,-self.inod,self.inod,40,-2.,2.)
        self.root.h2d(self.prefix+'aci',60,-self.inod,self.inod,40,-2.,2.)
        self.root.h2d(self.prefix+'fchi',60,-self.inod,self.inod,40,0.,10.)
        self.root.h2d(self.prefix+'rchi',60,-self.inod,self.inod,40,0.,10.)
        self.root.h2d(self.prefix+'chi',60,-self.inod,self.inod,40,0.,10.)
        self.root.h2d(self.prefix+'fchib',60,-self.inod,self.inod,40,0.,10.)
        self.root.h2d(self.prefix+'rchib',60,-self.inod,self.inod,40,0.,10.)
        self.root.h2d(self.prefix+'chib',60,-self.inod,self.inod,40,0.,10.)


        self.root.hprf(self.prefix+'faci_pf',60,-self.inod,self.inod,-2.,2.)
        self.root.hprf(self.prefix+'raci_pf',60,-self.inod,self.inod,-2.,2.)
        self.root.hprf(self.prefix+'aci_pf',60,-self.inod,self.inod,-2.,2.)
        self.root.hprf(self.prefix+'fchi_pf',60,-self.inod,self.inod,0.,10.)
        self.root.hprf(self.prefix+'rchi_pf',60,-self.inod,self.inod,0.,10.)
        self.root.hprf(self.prefix+'chi_pf',60,-self.inod,self.inod,0.,10.)
        self.root.hprf(self.prefix+'fchib_pf',60,-self.inod,self.inod,0.,10.)
        self.root.hprf(self.prefix+'rchib_pf',60,-self.inod,self.inod,0.,10.)
        self.root.hprf(self.prefix+'chib_pf',60,-self.inod,self.inod,0.,10.)
        return True

    def execute(self):
        vals = self.evt['rec/kfs/scan']
        if (not vals): return False
        inodes = vals.keys()
        inodes.sort()
        print 'inodes ',inodes
        for inode in inodes:
            rkf,fkf = vals[inode]
            if (not rkf or not fkf): continue
            rac = achi(rkf.nodes)[0]
            fac = achi(fkf.nodes)[0]
            rchi = chindf(rkf.nodes)
            fchi = chindf(fkf.nodes)
            rchib = chiblocks(rkf.nodes)
            print ' inode rchib ',inode,rchib
            fchib = chiblocks(fkf.nodes)
            print ' inode fchib ',inode,fchib
            if (len(rchib)>0): rchib=rchib[0]
            else: rchib=0.
            if (len(fchib)>0): fchib=fchib[0]
            else: fchib=0.
            print ' inode chi2 ',inode,rchi,fchi
            print ' inode ac   ',inode,rac,fac
            print ' inode chib   ',inode,rchib,fchib
            self.root.fill(self.prefix+'faci',inode,fac)
            self.root.fill(self.prefix+'raci',inode,rac)
            self.root.fill(self.prefix+'aci',inode,fac+rac)
            self.root.fill(self.prefix+'fchi',inode,fchi)
            self.root.fill(self.prefix+'rchi',inode,rchi)
            self.root.fill(self.prefix+'chi',inode,0.5*(rchi+fchi))
            self.root.fill(self.prefix+'fchib',inode,fchib)
            self.root.fill(self.prefix+'rchib',inode,rchib)
            self.root.fill(self.prefix+'chib',inode,0.5*(rchib+fchib))

            self.root.fill(self.prefix+'faci_pf',inode,fac)
            self.root.fill(self.prefix+'raci_pf',inode,rac)
            self.root.fill(self.prefix+'aci_pf',inode,fac+rac)
            self.root.fill(self.prefix+'fchi_pf',inode,fchi)
            self.root.fill(self.prefix+'rchi_pf',inode,rchi)
            self.root.fill(self.prefix+'chi_pf',inode,0.5*(rchi+fchi))
            self.root.fill(self.prefix+'fchib_pf',inode,fchib)
            self.root.fill(self.prefix+'rchib_pf',inode,rchib)
            self.root.fill(self.prefix+'chib_pf',inode,0.5*(rchib+fchib))

        return True


class PullEventDisplay(IAlg):

    def define(self):
        self.prefix = 'ed'
        self.nhits = 100
        self.ievt = 0
        return

    def execute(self):
        kfs = self.evt['rec/kfs']
        if (not kfs): return False
        self.ievt+=1

        for kf in kfs:

            self.root.h1d(self.prefix+'_z',self.nhits,0,self.nhits)
            self.root.h1d(self.prefix+'_ee',self.nhits,0,self.nhits)
            self.root.h1d(self.prefix+'_fchi',self.nhits,0,self.nhits)
            self.root.h1d(self.prefix+'_schi',self.nhits,0,self.nhits)
            self.root.h1d(self.prefix+'_xpull',self.nhits,0,self.nhits)
            self.root.h1d(self.prefix+'_ypull',self.nhits,0,self.nhits)
            self.root.h1d(self.prefix+'_txpull',self.nhits,0,self.nhits)
            self.root.h1d(self.prefix+'_typull',self.nhits,0,self.nhits)
            
            nn = len(kf.nodes)
            for i in range(nn):
                node = kf.nodes[i]
                rx,ry,rtx,rty,rene = node.getstate('true').vec
                mtype = 'smooth'
                x,sx = node.param(mtype,0)
                y,sy = node.param(mtype,1)
                tx,stx = node.param(mtype,2)
                ty,sty = node.param(mtype,3)
                ee,see = node.param(mtype,4)
                self.root.fill(self.prefix+'_z',i,node.zrun)
                self.root.fill(self.prefix+'_e',i,rene)
                self.root.fill(self.prefix+'_fchi',i,node.getchi2('filter'))
                self.root.fill(self.prefix+'_schi',i,node.getchi2('smooth'))
                self.root.fill(self.prefix+'_xpull',i,(x-rx)/sx)
                self.root.fill(self.prefix+'_ypull',i,(y-ry)/sy)
                self.root.fill(self.prefix+'_txpull',i,(tx-rtx)/stx)
                self.root.fill(self.prefix+'_typull',i,(ty-rty)/sty)
            
            names = map(lambda x: self.prefix+'_'+x,['z','e','fchi','schi','xpull','ypull','txpull','typull'])
            hs = map(lambda name: self.root.get(name),names)

            t0 = tcanvas(hs[:4],nx=2,name='evtdis_chi2_evt'+str(self.ievt))
            t1 = tcanvas(hs[4:],nx=2,name='evtdis_pull_evt'+str(self.ievt))
            raw_input('press key ')    

        return True


#-----------------------------
# Alex
#-----------------------------

def dovertex():

    alex = Alex()
    alex.nevts = 1
    root = ROOTSvc('root','akfbeta.root')
    alex.addsvc(root)

    # simulate states
    agenstates = GenerateDoubleBeta('agenstates')
    #agenstates = GenerateBeta('agenstates')
    alex.addalg(agenstates)

    # histogram simulated states
    hisgenstates = HistosGenerateStates('hisgenstates')  
    hisgenstates.imports.append('root')
    alex.addalg(hisgenstates)

    agendigits = CreateDigits('agendigits')
    agendigits.dz = 0.5
    alex.addalg(agendigits)  

    #agennodes = CreateNodes('agennodes')
    #alex.addalg(agennodes)

    betanodes = DoubleBetaNodes('betanodes')
    alex.addalg(betanodes)

    akf = FitScanNodes('ascan')
    #akf.path = 'rec/nodes'
    #akf.opath = 'rec/kfs'
    alex.addalg(akf)

    hakf = HistosScanNodes('hscan')
    hakf.imports.append('root')
    alex.addalg(hakf)

    alex.run()
    
def do():

    alex = Alex()
    alex.nevts = 1
    root = ROOTSvc('root','akfbeta.root')
    alex.addsvc(root)

    # simulate states
    agenstates = GenerateDoubleBeta('agenstates')
    #agenstates = GenerateBeta('agenstates')
    alex.addalg(agenstates)

    # histogram simulated states
    hisgenstates = HistosGenerateStates('hisgenstates')  
    hisgenstates.imports.append('root')
    alex.addalg(hisgenstates)

    agendigits = CreateDigits('agendigits')
    agendigits.dz = 0.5
    alex.addalg(agendigits)  

    #agennodes = CreateNodes('agennodes')
    #alex.addalg(agennodes)

    betanodes = DoubleBetaNodes('betanodes')
    alex.addalg(betanodes)

    akf = FitNodes('akf')
    akf.path = 'rec/nodes'
    akf.opath = 'rec/kfs'
    alex.addalg(akf)

    hisakf = HistosKFFit('hisakf')
    hisakf.imports.append('root')
    hisakf.prefix='kf_'
    hisakf.path='rec/kfs'
    alex.addalg(hisakf)

    arevnodes = ReverseNodes('arevnodes')
    arevnodes.path = 'rec/nodes'
    arevnodes.opath = 'rec/nodes/rev'
    alex.addalg(arevnodes)

    akfrev = FitNodes('akfrev')
    akfrev.path = 'rec/nodes/rev'
    akfrev.opath = 'rec/kfs/rev'
    alex.addalg(akfrev)

    hisakfrev = HistosKFFit('hisakfrev')
    hisakfrev.imports.append('root')
    hisakfrev.prefix='kfrev_'
    hisakfrev.path='rec/kfs/rev'
    alex.addalg(hisakfrev)

#    anarev = AnaReverse('anrev')
#    anarev.imports.append('root')
#    alex.addalg(anarev)

    alex.run()    
    
if __name__ == '__main__':

    dovertex()
