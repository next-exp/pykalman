"""
functions and algorithms to generate using a toy MC and fit with the Kalman Filter.
single electrons and double-beta events.

The program run using the alex framework

"""


from math import *
import random
from copy import deepcopy
from array import array

from alex import IAlg
from alex import Alex
from rootsvc import ROOTSvc,ttree

from KFBase import KFVector, KFMatrix, KFMatrixNull, KFMatrixUnitary
from kfnext import NEXT, nextgenerator, nextfilter, V0, H0, simplegenerator, simplefilter
from kfgenerator import zavesample,zrunsample
from kffilter import randomnode, KFData, KFNode
from kfzline import zstate
from messenger import Message
from ROOT import TCanvas,TPolyLine3D,TTree,TGraph
from KFTrack import *

#-------------------------------
# Functions
#-------------------------------


# set the condition of NEXT
gpgas = 15.
gnext = NEXT(gpgas)
gnextgen = nextgenerator(gnext)
gnextfit = nextfilter()

#radlen = 1500.
#gnextgen = simplegenerator(radlen)
#gnextfit = simplefilter(radlen)

MTYPE = 'filter'
NBLOCK = 5
#MTYPE = 'smooth'

def genele(ene = 2.5):
    """ generate an electron with a kinection energy (ene).
    return a list of states (x,y,z,ux,uy,ux,ene), 
    where ux,uy,uz are the cos-director and ene is the kinetic energy
    """
    theta = acos(random.uniform(-1.,1.))
    phi = random.uniform(0,2.*pi)
    ux,uy,uz = sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)
    #ux,uy,uz=0.,0.,1.
    state0 = (0.,0.,0.,ux,uy,uz,ene)
    states = gnextgen.generate(state0)
    #print 'genele-states-'states
    return states

def genbeta(ene = 2.5, emin=0.20):
    """ generate a double beta (two electrons) with total kinetic energy (ene).
    Return a list with two list of states with (x,y,z,ux,uy,uz,ene)
    where ux,uy,uz are the cos-director and ene is the kinetic energy
    """
    ei0 = random.uniform(emin,ene-emin)
    ei1 = ene-ei0
    e0 = min(ei0,ei1)
    e1 = max(ei0,ei1)
    states1 = genele(e0)
    states2 = genele(e1)
    #print 'genbeta-states1 ',states1
    #print 'genbeta-states2 ',states2
    
    return [states1,states2]

def createzdigits(states,zs):
    """ from the list os states, create a list of digits.
    an state is a list with (x,y,y,ux,uy,uz,ene)
    a digit has (x,y,z,delta-ene). 
    return the list of digits and the states used to generate the digits
    """
    zstates = zavesample(states)
    denes = map(lambda i: zstates[i-1][-1]-zstates[i][-1],range(1,len(zstates)))    
    denes[-1]= denes[-1]+(zstates[-1][-1]) # add the last deposition of energy to the last state
    #print "createzdigits denes "
    #print "createzdigits sum(denes) ",sum(denes)
    digits = map(lambda st,dene: [st[0],st[1],st[2],dene],zstates[1:],denes)
    #print 'createzdigits digits'
    #note that it removes the first state!
    return (digits,zstates[1:])

def createhits(digits,xres=0.1):  #Resolution Change
    """ from the list of digits (x,y,z,delta-ene), return a list of hits.
    Each hit is computed with a given resolution
    A hit is a KFData with a vector (x,y) and a cov-matrix. 
    Hit also has an attribute with the delta-ene (dene)
    """
    hits = []
    for digit in digits:
        x,y,z,dene = digit
        x = x+random.gauss(0.,xres)
        y = y+random.gauss(0.,xres)
        V = (xres*xres)*KFMatrixUnitary(2)
        hit = KFData(KFVector([x,y]),V,zrun=z)
        hit.dene = dene
        hits.append(hit)
    #print 'createhits ',hits
    return hits

def createnodes(hits,H=H0,states=None):
    """ from the list of hits (a KFData with (x,y) and cov-matrix), 
    generate a list of nodes using the proyection matrix H0.
    if states are provided they are stored as 'true' states
    """
    nodes = map(lambda hit: KFNode(hit,H),hits)
    if (states):
        zstates = map(zstate,states)
        for i,node in enumerate(nodes):
            #print zstates[i]
            node.setstate('true',zstates[i])
    #print "create nodes "
    return nodes

def preparenodes(nodes):
    """ this method is called before the kalman-filter to prepare the nodes.
    It add the delta-energy of the nodes in the order of propagation, to used during the fit.
    That is, the first node has the total energy of the track.
    """
    denes = map(lambda nd: nd.hit.dene,nodes)
    for i,node in enumerate(nodes):
        node.hit.ene = sum(denes[i:])
    enes = map(lambda nd: nd.hit.ene,nodes)
    zs = map(lambda nd: nd.zrun,nodes)
    #print ' preparenodes enes ',enes
    #print ' zs ',zs
    return

def seedstate(nodes):
    """ from a list of nodes generate the seed state (a KFData (x,y,tx,ty,ene) and cov-matrix)
    where tx,ty are the tangent in the x,y axis
    """
    x0 = nodes[0].hit.vec
    x1 = nodes[1].hit.vec
    z0 = nodes[0].zrun
    dz = nodes[1].zrun-z0
    if (dz==0.): 
        print ' Fit failed, not valid input nodes ',z0,dz
        #print ' node 0 ',nodes[0]
        #print ' node 1 ',nodes[1]
        return False,None
    ene = nodes[0].hit.ene
    tx,ty = (x0[0]-x1[0])/dz,(x0[1]-x1[1])/dz
    xx = KFVector([x0[0],x0[1],tx,ty,ene])
    zdir = dz/abs(dz)
    C1 = 100.*KFMatrixUnitary(5)
    state0 = KFData(xx,C1,z0-0.1*zdir,pars={'uz':zdir})
    #print "seedstate ",state0
    return state0

def fitnodes(nodes):
    """ fit the nodes using the kalman filter, 
    it return the result of the fit: (cc,kk):
        cc is list with (chi2-filter, chi2-smooth), 
        kf is the kalman filter object with the full information 
    """
    #print ' fitnodes ',len(nodes)
    preparenodes(nodes)
    state0 = seedstate(nodes)
    if state0 == (False,None): return (False,False,False),False
    gnextfit.clear()
    gnextfit.setnodes(nodes)
    cc = gnextfit.fit(state0)
    kf = deepcopy(gnextfit)
    #print 'fitnodes ',cc    
    return cc,kf

def chiblocks(nodes,maxchi=30.,mtype=MTYPE,nblock=NBLOCK):
    """ return the chi2 of the track divided in blocks of (nblock) nodes, 
    removing the nodes with large chi2 (maxchi),
    mtype can be 'filter' of 'smooth'
    """
    if (mtype=='filter'): nodes = nodes[2:]
    chis = map(lambda nd: nd.getchi2(mtype),nodes)
    chis = filter(lambda chi: chi<maxchi,chis)
    nn = len(chis)
    nb = int(nn/(1.*nblock))
    bchis = map(lambda i: sum(chis[i*nblock:(i+1)*nblock])/(2.*nblock-4.),range(nb))
    if (len(bchis)<=0): bchis=[0.,]
    #print 'chiblocks ',nn,nb,bchis
    #chis.append(sum(chi2[nblock*nb:])/(2*nblock-4))
    return bchis

def chindf(nodes,maxchi=30.,mtype=MTYPE):
    """ computes the chi/ndf of the nodes, removing nodes with large chi2 (maschi2),
    mtype can be 'filter' or 'smooth'
    """
    chi = map(lambda nd: nd.getchi2(mtype),nodes)
    chi = filter(lambda ch: ch<maxchi,chi)
    nn = len(chi)
    if nn<=2: chi = 1.
    else: chi = sum(chi)/(2.*nn-4.)
    return chi

def achi(nodes,maxchi=30.,mtype=MTYPE,nblock=NBLOCK):
    """ return the chi2-asymmetry.
    That is the asymmetry of the first and last block of the track.
    the number of nodes in each block is (nblock).
    nodes with large chi2 (maxchi) are removed.
    mtype can be 'filter' or 'smooth'
    """
    chis = chiblocks(nodes,maxchi,mtype,nblock)
    if (len(chis)<=1): return 0.
    chi1,chi2 = chis[0],chis[-1]
    ac = (chi1-chi2)/(chi1+chi2)
    schi = sum(chis)/len(chis)
    c1 = chi1/schi
    #print ' achi ',ac,c1
    return ac

def achiv0(nodes,maxchi=30.,m=3):
    """ return the chi2-asymmetry.
    the nodes are divided in m-parts and the assymmetry is computed from the first and last part.
    nodes with large chi2 (maxchi) are removed.
    """
    nn = len(nodes)
    ni = int(nn/m)
    ac = 0.
    if (ni<=2): return ac
    for mtype in ['filter','smooth']:
        n0 = 0
        if (mtype == 'filter'): n0 = 2
        ch1 = chindf(nodes[n0:ni],mtype=mtype)
        ch3 = chindf(nodes[(m-1)*ni:],mtype=mtype)
        ac = (ch1-ch3)/(ch1+ch3)
    # print "achi ",acs
    return ac

def foms(nodes):
    """ return figure of merit of the fit (fom)
    (chi,chib,ac): total chi/ndf, chi2 of the first block, and the chi2-asymmetry
    """
    chi = chindf(nodes)
    chib = chiblocks(nodes)[0]
    ac = achiv0(nodes)
    return chi,chib,ac

#------------------------
#    algorithms
#------------------------

class GenerateBeta(IAlg):
    """ Algorithm to generate states (x,y,z,ux,uy,uz,p)
    """
    def define(self):
        self.E0 = 2.5 # MeV
        return
        
    def initialize(self):
        """ book the ntuple and put it in ROOTSvc
        """
        n = 492
        labels = [('x','F',n),('y','F',n),('z','F',n),
                  ('ux','F',n),('uy','F',n),('uz','F',n),
                  ('ee','F',n)]
        tree = ttree(self.name,labels)
        self.root.put(tree)
        return True

    def execute(self):
        states = [genele(self.E0),]
        track = states[0]

        x = map(lambda s: s[0],track)
        y = map(lambda s: s[1],track)
        z = map(lambda s: s[2],track)
        ux = map(lambda s: s[3],track)
        uy = map(lambda s: s[4],track)
        uz = map(lambda s: s[5],track)
        ee = map(lambda s: s[6],track)
        tup = {'x':x,'y':y,'z':z,'ux':ux,'uy':uy,'uz':uz,'ee':ee}
        self.root.fill(self.name,tup)

        
        self.evt['sim/states'] = states
        val = map(lambda st: (len(st),st[0],st[-1]),states)
        self.msg.verbose(self.name,val)        
        return True

class GenerateDoubleBeta(IAlg):
    """ Algorithm to generate a doble beta event
    """

    def define(self):
        self.E0 = 2.5 # MeV
        return

    def execute(self):
        states = genbeta(self.E0)
        self.evt['sim/states'] = states
        data = map(lambda seg: (len(seg),seg[0],seg[-1]),states)
        self.msg.verbose(self.name,' states ',data)
        return True

class CreateDigits(IAlg):
    """ Algorithm to generate digits, separated a distance in z (dz)
    """

    def define(self):
        self.dz = 0.5
        self.ZZ = 30.
        return

    def initialize(self):
        n = int(2.*self.ZZ/self.dz)
        self.zs = map(lambda i: -1.*self.ZZ+i*self.dz,range(n))
        self.msg.verbose(self.name,' zs ',self.zs)
        return True

    def execute(self):
        states = self.evt['sim/states']
        if (not states): return False
        digs,zsts= [],[]
        for seg in states:
            dig,zst = createzdigits(seg,self.zs)
            digs.append(dig)
            zsts.append(zst)
        ok = len(digs)>0
        #print digs
        #print zsts
        if (ok): 
            self.evt['sim/digits']=digs
            self.evt['sim/zstates']=zsts
        data1 = map(lambda dig: (len(dig),dig[0],dig[-1]),digs)
        data2 = map(lambda zst: (len(zst),zst[0],zst[-1]),zsts)
        self.msg.verbose(self.name,' digits ',data1)
        self.msg.verbose(self.name,' zstats ',data2)    
        return True


class CreateNodes(IAlg):
    """ Algorithm to create hits and nodes.
    The resolution of the measurement is controlled by xres
    """

    def define(self):
        self.xres = 0.1             #Resolution Change
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
        self.msg.verbose(self.name,' hits ',dat1)
        #self.msg.verbose(self.name,'nodes ',dat2)
        return ok

class ReverseNodes(IAlg):
    """ Algorithm to reverse the nodes
    """

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
        self.msg.verbose(self.name,' nodes ',val)
        return True

class DoubleBetaNodes(IAlg):
    """ Algorithm to create the nodes of one track from a doble-beta event!
    """

    def execute(self):
        digs = self.evt['sim/digits'][0]
        zsts = self.evt['sim/zstates'][0]
#        print 'digs',len(digs),len(digs[0]),len(digs[1])#digs
        if (not digs or not zsts): return False
        #print 'digs ',len(digs),len(digs[0]),len(digs[1])
        if (len(digs)!=2): return False
        dig = deepcopy(digs[0][0][1:])
        zst = deepcopy(zsts[0][0][1:])
        dig.reverse(); zst.reverse()
        dig+=deepcopy(digs[1][0]); zst+=deepcopy(zsts[1][0])
        #print ' digs - total ',len(dig)
#        print 'd',len(dig),dig
        hits = [createhits(dig),]
        nods = [createnodes(hits[0],H0,zst),]
        ok = len(nods)>0
#        print 'ok',ok
        if (ok):
            self.evt['rec/hits'] = hits
            self.evt['rec/nodes'] = nods
        dat1 = map(lambda hit: (len(hit),hit[0]),hits)
        #dat2 = map(lambda nod: (len(nod),nod[0]),nods)
        self.msg.verbose(self.name,'hits ',dat1)
        #self.msg.verbose(self.name,'nodes ',dat2)
        return ok

class FitNodes(IAlg):
    """ Algorithm to fit the nodes
    """

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
            self.msg.verbose(self.name,' fit ',cc)
            if ok: kfs.append(kf)
        ok = len(kfs)>0
        if (ok): self.evt[self.opath]=kfs
        vals = map(lambda kf: (kf.cleanchi2('filter'),kf.cleanchi2('smooth')),kfs)
        self.msg.info(self.name,' chi2 ',vals)
        return ok

class FitNodes2(IAlg):
    """ Algorithm to fit the nodes
    """
    def Invert(self,lis):
        m = max(lis)
        return [m-i for i in lis]
        
    def Plot4D( self, x, y, z, t, markerstyle = 20, markersize = 1 ):
        '''
            Plot a 3D dataset (x,y,z) with an extra color coordinate (t).
        '''
        data = array( 'd', [0.] * 4 )
        tree = TTree('DummyTree','DummyTree')
        tree.Branch('xyzt', data, 'x/D:y:z:t')

        for datai in zip(x,y,z,t):
            data[0], data[1], data[2], data[3] = datai
            tree.Fill()
        tree.SetMarkerStyle( markerstyle )
        tree.SetMarkerSize( markersize )
        c = TCanvas()
        tree.Draw('x:y:z:t','','zcol')
        return c, tree
        
    def PlotTrack( self, xt, yt, zt, Et, kf, state = 'filter' ):
        
        xf = [ node.getstate( state ).vec[0] for node in kf.nodes ]
        yf = [ node.getstate( state ).vec[1] for node in kf.nodes ]
        zf = [ node.zrun for node in kf.nodes ]
        Ef = self.Invert( [ node.getstate( state ).vec[-1] for node in kf.nodes ] )
        chi2 = [ node.chi2[state] for node in kf.nodes ]
        line = TPolyLine3D( len(xf), array('f',zf), array('f',yf), array('f',xf) )
        a = self.Plot4D( xt, yt, zt, Et )
        line.Draw('same')
    #    cc = TCanvas()
    #    g = Graph( range(len(xf)), chi2, markerstyle = 20 )
    #    g.Draw('AP')
        g, cc = 0, 0
        return a, line,g,cc

    def define(self):
        self.path = 'rec/nodes'
        self.opath = 'rec/kfs'
        return True

    def execute(self):
        nodes = self.evt[self.path][0]
#        print len(nodes),nodes[0]
        if (not nodes): return False
        kfs = []
        cc,kf = fitnodes(nodes)
        ok,fchi,schi = cc
        self.msg.verbose(self.name,' fit ',cc)
        #print 'ok1',ok
        if ok: kfs.append(kf)
        ok = len(kfs)>0
        #print 'ok2',ok
        if (ok): self.evt[self.opath]=kfs
        vals = map(lambda kf: (kf.cleanchi2('filter'),kf.cleanchi2('smooth')),kfs)
        self.msg.info(self.name,' chi2 ',vals)
#       print kfs
        #plot 
        if not ok:
            return ok
            
        # Descomentando todo lo siguiente tenemos un plot 3d para cada traza con el ajuste del KF
        #states = self.evt['sim/zstates'][0]
        #x,y,z,ux,uy,uz,E = zip(*states)
        #print 'kfs',type(kf),len(kf.nodes)
        #l=self.PlotTrack(x,y,z,E,kf)
        #raw_input('enter')
        
        return ok

class FitScanNodes(IAlg):
    """ Algorithm to do a scan along the track and fit in both directions
    """

    def define(self):
        self.path = 'rec/nodes/bb'
        self.opath = 'rec/kfs/scan'
        self.maxscan = 40
        self.iscan = 1 # step in the scan (1=each node)
        self.nmin = 2 # minimun number of nodes to do a fit
        return

    def execute(self):
        segs = self.evt[self.path]
        digs = self.evt['sim/digits']
        if (not segs or not digs): return False
        nver = len(digs[0])
        nodes = segs[0]
        ntot = len(nodes)-2.*self.nmin
        inodes = map(lambda i: self.nmin+self.iscan*i,range(1,int(ntot/self.iscan)-1))
        if (not (nver in inodes)): inodes.append(nver)
        inodes.sort()
        inodes = filter(lambda i: abs(i-nver)<self.maxscan,inodes)
        #print 'nnodes nvertex',ntot,nver
        #print 'inodes ',inodes
        vals = {}
        for inode in inodes:
            rnodes = deepcopy(nodes[:inode])
            rnodes.reverse()
            fnodes = deepcopy(nodes[inode:])
            #print ' forward ',inode
            fcc,fkf = fitnodes(fnodes)
            #print ' reverse ',inode
            rcc,rkf = fitnodes(rnodes)
            vals[inode-nver]=(rkf,fkf)
            self.msg.verbose('scan i ',inode-nver)
            self.msg.verbose('scan r ',rcc,achi(rkf.nodes))
            self.msg.verbose('scan f ',fcc,achi(fkf.nodes))
        self.evt['rec/kfs/scan']=vals
        return True

#-----------------------------------
#   histogramming!!
#------------------------------------

class HistosGenerateStates(IAlg):
    """ Algorithm to histogram the generated states
    """

    def define(self):
        self.nstates = 160 # max number of states
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
    """ Algorithm to histogram the kalfim filter results
    """

    def define(self):
        self.nhits = 160 # max number of hits on track
        self.maxchi2 = 30.
        self.rev = False
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
        self.root.h1d(self.prefix+'fchieeLE',100,0.,10.) #######################
        self.root.h1d(self.prefix+'fchieeHE',100,0.,10.) #######################
#        self.root.h1d(self.prefix+'nhits',100,-5.,5.)

        self.root.h2d(self.prefix+'xi',self.nhits,0,self.nhits,100,-10.,10.)
        self.root.h2d(self.prefix+'yi',self.nhits,0,self.nhits,100,-10.,10.)
        self.root.h2d(self.prefix+'zi',self.nhits,0,self.nhits,100,-5.,20.)
        self.root.h2d(self.prefix+'txi',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'tyi',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'eei',self.nhits,0,self.nhits,100,0.,3.)
        self.root.h2d(self.prefix+'fchii',self.nhits,0,1.,100,0.,10.)
        self.root.h2d(self.prefix+'fchiee',100,0,2.5,100,0.,10.) #######################
        self.root.h2d(self.prefix+'schii',self.nhits,0,1.,100,0.,10.)
        self.root.hprf(self.prefix+'fchii_pf',self.nhits,0,1,0.,self.nhits)
        self.root.hprf(self.prefix+'schii_pf',self.nhits,0,1,0.,self.nhits)
        self.root.h2d(self.prefix+'xpulli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'ypulli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'txpulli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'typulli',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d(self.prefix+'fchiHELE',100,0.,10.,100,0.,10.)
    def execute(self):
        kfs = self.evt[self.path]
        #print kfs
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
        fac = achi(kf.nodes,mtype='filter')
        sac = achi(kf.nodes,mtype='smooth')
        self.root.fill(self.prefix+'fac',fac)
        self.root.fill(self.prefix+'sac',sac)
        
    	eaux = 2.5/2.
        if self.rev:
			eaux = 2.5-eaux
        chiauxHE = 0.
        chiauxLE = 0.
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

            root.fill(self.prefix+"xi",float(i)/nhits,x)
            root.fill(self.prefix+"yi",float(i)/nhits,y)
            root.fill(self.prefix+"zi",float(i)/nhits,rz)
            root.fill(self.prefix+"txi",float(i)/nhits,tx)
            root.fill(self.prefix+"tyi",float(i)/nhits,ty)
            root.fill(self.prefix+"eei",float(i)/nhits,ee)
            
            
            
            
            if (schi<=self.maxchi2):
                root.fill(self.prefix+'xpulli',float(i)/nhits,(x-rx)/sx)
                root.fill(self.prefix+'ypulli',float(i)/nhits,(y-ry)/sy)
                root.fill(self.prefix+'txpulli',float(i)/nhits,(tx-rtx)/stx)
                root.fill(self.prefix+'typulli',float(i)/nhits,(ty-rty)/sty)
                root.fill(self.prefix+'schii',float(i)/nhits,schi)
                root.fill(self.prefix+'schii_pf',float(i)/nhits,schi)
            if (fchi<=self.maxchi2):
                root.fill(self.prefix+'fchii',float(i)/nhits,fchi)
                root.fill(self.prefix+'fchii_pf',float(i)/nhits,fchi)
                root.fill(self.prefix+'fchiee',ee,fchi) ###
            
                if ee >= eaux:
                    chiauxHE += fchi
                    aux = i
                else:
                
                    chiauxLE += fchi
            
            #print ' kf-node i, z, ene, fchi ',i,rz,ree,fchi
        #print chiauxHE/aux
        #print chiauxLE/aux
        root.fill(self.prefix+'fchieeLE',chiauxLE)
        root.fill(self.prefix+'fchieeHE',chiauxHE)
        return True


class HistosAnaReverse(IAlg):
    """ Algorithm to histogram the fits in the forward and reverse mode
    """

    def define(self):
        self.nhits = 160 # max number of hits on track
        self.maxchi2 = 30.
        self.prefix = 'ana_'
        self.pathnor = 'rec/kfs'
        self.pathrev = 'rec/kfs/rev'
        self.graph = TGraph()
        return

    def initialize(self):
        self.root.h2d(self.prefix+'chi',20,0.,5.,20,0.,5.)
        self.root.h2d(self.prefix+'chib',20,0.,12.,20,0.,12.)
        self.root.h2d(self.prefix+'ac',20,-1.,1.,20,-1.,1.)
        self.root.h3d(self.prefix+'xiFB',1000,0.,30,1000,0.,30.,100,0,self.nhits)    
        self.root.h2d(self.prefix+'meanChi',100,0.,10.,100,0.,10.)  
        self.root.h2d(self.prefix+'AsimChii',100,0.,1,100,-1,1)   
        self.root.h2d(self.prefix+'halfAsimChi',100,-1.,1,100,-1,1)  
        self.root.h2d(self.prefix+'halfMeanIn',100,0.,10,100,0,10)  
        self.root.h1d(self.prefix+'AsimIn',100,-1.,1.)  
        self.root.h1d(self.prefix+'AsimMed',100,-1.,1.) 
        self.root.h1d(self.prefix+'AsimFin',100,-1.,1.)  
        self.root.h1d(self.prefix+'AsimAsim',100,-1.05,1.05) 
        self.root.h2d(self.prefix+'AsimBoth',100,-1.,1.,100,-1.,1.) 
        return True

    def execute(self):
        
        nkfs = self.evt[self.pathnor]
        rkfs = self.evt[self.pathrev]
        
        if (not nkfs or not rkfs): return False
        for nkf in nkfs:
            for rkf in rkfs:
                    nhits = min(len(rkf),len(nkf))
                    
                    auxn = 0.
                    auxr = 0.
                    norIn = 0.
                    norFin = 0.
                    norMed = 0.
                    revIn = 0.
                    revFin = 0.
                    revMed = 0.
                    
                    asimI = 0.
                    asimM = 0.
                    asimF = 0.
  
                    
                    for i in range(nhits):
                        paso = float(i)/nhits
                        nnode = nkf.nodes[i]
                        nfchi = nnode.getchi2('filter')
                        auxn += nfchi
                        nschi = nnode.getchi2('smooth')
                        rnode = rkf.nodes[i]
                        rfchi = rnode.getchi2('filter')
                        auxr += rfchi
                        rschi = rnode.getchi2('smooth')
                        self.root.fill(self.prefix+'xiFB',nfchi,rfchi,i)
                        asim = (nfchi-rfchi)/(nfchi+rfchi)
                        self.root.fill(self.prefix+'AsimChii',paso,asim)
                        #self.graph.SetPoint(i,nfchi,rfchi)
                        if paso <= 0.33:
                            norIn += nfchi
                    	    revIn += rfchi
                    	    asimI += asim
                    	    ni = i
                    	elif paso >= 0.66:
                    	    norFin += nfchi
                    	    revFin += rfchi
                    	    asimF +=asim
                          
                    	else:
                    	    norMed += nfchi
                    	    revMed += rfchi
                    	    asimM += asim
                    	    nm = i
                    
                    nm -= ni	    
                    self.root.fill(self.prefix+'meanChi',auxn/nhits,auxr/nhits)
                    norIn *= 1./ni
                    revIn *= 1./ni
                    asimI *= 1./ni
                    norMed *= 1./nm
                    revMed *= 1./nm
                    asimM *= 1./nm
                    norFin *= 1./(nhits-ni-nm)
                    revFin *= 1./(nhits-ni-nm)
                    asimF *= 1./(nhits-ni-nm)
                    
                    asimNorm = (norIn-norMed)/(norIn+norMed)
                    asimRev = (revIn-revMed)/(revIn+revMed)
                    asimBIn = (revIn-norIn)/(revIn+norIn)
                    asimBMed = (revMed-norMed)/(revMed+norMed)
                    asimBFin = (revFin-norFin)/(revFin+norFin)
                    asimasim = (-asimI+asimF)/(abs(asimI)+abs(asimF))
                    
                    self.root.fill(self.prefix+'halfAsimChi',asimNorm,asimRev)
                    self.root.fill(self.prefix+'halfMeanIn',norIn,revIn)
                    
                    self.root.fill(self.prefix+'AsimIn',asimI)
                    self.root.fill(self.prefix+'AsimMed',asimM)
                    self.root.fill(self.prefix+'AsimFin',asimF)
                    self.root.fill(self.prefix+'AsimAsim',asimasim)
                    self.root.fill(self.prefix+'AsimBoth',asimBIn,asimBFin)
                    fi = open('Single_Irene_15_0T_10.dat','a')
                    liste = [auxn/nhits, auxr/nhits, norIn, revIn, asimI,norMed,revMed, asimM,norFin,revFin, asimF, asimasim, asimNorm, asimRev, asimBIn, asimBFin,asimBMed]
                    fi.write(' '.join(map(str,liste))+'\n')
                    fi.close()
        pars = zip(nkfs,rkfs)
        i = 0
        for par in pars:
            
            rk,fk = par
            fchi,fchib,fac = foms(fk.nodes)
            rchi,rchib,rac = foms(rk.nodes)
            self.root.fill(self.prefix+'chi',fchi,rchi)
            self.root.fill(self.prefix+'chib',fchib,rchib)
            self.root.fill(self.prefix+'ac',fac,rac)
            
            
            i += 1  
        return True
        
    def finalize(self):
        #c = TCanvas()
        x=1
        #self.graph.Draw('AP')
        #c.SaveAs('graph.root')
            

class HistosScanNodes(IAlg):
    """ Algorithm to histogram the scan of the track
    """

    def define(self):
        self.prefix = 'sn_'
        self.path = 'rec/kfs/scan'
        self.inod = 10
        return

    def initialize(self):

        self.root.h1d(self.prefix+'vtx_fchi',50,0.,5.)
        self.root.h1d(self.prefix+'vtx_rchi',50,0.,5.)
        self.root.h1d(self.prefix+'vtx_fchib',50,0.,10.)
        self.root.h1d(self.prefix+'vtx_rchib',50,0.,10.)
        self.root.h1d(self.prefix+'vtx_fac',50,-1.,1.)
        self.root.h1d(self.prefix+'vtx_rac',50,-1.,1.)
        self.root.hprf(self.prefix+'vtx_fchii_pf',30,0.,30,0.,10.)
        self.root.hprf(self.prefix+'vtx_rchii_pf',30,0.,30,0.,10.)

        self.root.h1d(self.prefix+'nvtx_fchi',50,0.,5.)
        self.root.h1d(self.prefix+'nvtx_rchi',50,0.,5.)
        self.root.h1d(self.prefix+'nvtx_fchib',50,0.,10.)
        self.root.h1d(self.prefix+'nvtx_rchib',50,0.,10.)
        self.root.h1d(self.prefix+'nvtx_fac',50,-1.,1.)
        self.root.h1d(self.prefix+'nvtx_rac',50,-1.,1.)
        self.root.hprf(self.prefix+'nvtx_fchii_pf',40,0.,40,0.,10.)
        self.root.hprf(self.prefix+'nvtx_rchii_pf',40,0.,40,0.,10.)

        self.root.h2d(self.prefix+'vnv_fchi',20,0.,5.,20,0.,5.)
        self.root.h2d(self.prefix+'vnv_fchib',20,0.,12.,20,0.,12.)
        self.root.h2d(self.prefix+'vnv_fac',20,-1.,1.,20,-1.,1.)
        self.root.h2d(self.prefix+'vnv_rchi',20,0.,5.,20,0.,5.)
        self.root.h2d(self.prefix+'vnv_rchib',20,0.,12.,20,0.,12.)
        self.root.h2d(self.prefix+'vnv_rac',20,-1.,1.,20,-1.,1.)

        self.root.h2d(self.prefix+'faci',60,-self.inod,self.inod,40,-2.,2.)
        self.root.h2d(self.prefix+'raci',60,-self.inod,self.inod,40,-2.,2.)
        self.root.h2d(self.prefix+'aci',60,-self.inod,self.inod,40,-2.,2.)
        self.root.h2d(self.prefix+'fchii',60,-self.inod,self.inod,40,0.,10.)
        self.root.h2d(self.prefix+'rchii',60,-self.inod,self.inod,40,0.,10.)
        self.root.h2d(self.prefix+'fchibi',60,-self.inod,self.inod,40,0.,10.)
        self.root.h2d(self.prefix+'rchibi',60,-self.inod,self.inod,40,0.,10.)
        self.root.h2d(self.prefix+'chibi',60,-self.inod,self.inod,40,0.,10.)

        self.root.hprf(self.prefix+'fac_pf',60,-self.inod,self.inod,-2.,2.)
        self.root.hprf(self.prefix+'rac_pf',60,-self.inod,self.inod,-2.,2.)
        self.root.hprf(self.prefix+'ac_pf',60,-self.inod,self.inod,-2.,2.)
        self.root.hprf(self.prefix+'fchi_pf',60,-self.inod,self.inod,0.,10.)
        self.root.hprf(self.prefix+'rchi_pf',60,-self.inod,self.inod,0.,10.)
        self.root.hprf(self.prefix+'chi_pf',60,-self.inod,self.inod,0.,10.)
        self.root.hprf(self.prefix+'fchib_pf',60,-self.inod,self.inod,0.,10.)
        self.root.hprf(self.prefix+'rchib_pf',60,-self.inod,self.inod,0.,10.)
        self.root.hprf(self.prefix+'chib_pf',60,-self.inod,self.inod,0.,10.)

        self.root.h2d(self.prefix+'v0a_fchii',60,-self.inod,self.inod,40,-1.,1.)
        self.root.h2d(self.prefix+'v0a_rchii',60,-self.inod,self.inod,40,-1.,1.)
        self.root.h2d(self.prefix+'v0a_fchibi',60,-self.inod,self.inod,40,-1.,1.)
        self.root.h2d(self.prefix+'v0a_rchibi',60,-self.inod,self.inod,40,-1.,1.)

        self.root.hprf(self.prefix+'v0a_fchi_pf',60,-self.inod,self.inod,-1.,1.)
        self.root.hprf(self.prefix+'v0a_rchi_pf',60,-self.inod,self.inod,-1.,1.)
        self.root.hprf(self.prefix+'v0a_fchib_pf',60,-self.inod,self.inod,-1.,1.)
        self.root.hprf(self.prefix+'v0a_rchib_pf',60,-self.inod,self.inod,-1.,1.)

        return True

    def execute(self):
        vals = self.evt['rec/kfs/scan']
        if (not vals): return False
        inodes = vals.keys()
        inodes.sort()
        #print 'inodes ',inodes

        def ifoms(inode):
            rkf,fkf = vals[inode]
            rvals = [0.,0.,0.]
            if (rkf): rvals = foms(rkf.nodes)
            fvals = [0.,0.,0.]
            if (fkf): fvals = foms(fkf.nodes)
            return (rvals,fvals)

        vs = ifoms(0)
        vrchi,vrchib,vrac = vs[0]
        vfchi,vfchib,vfac = vs[1]

        for inode in inodes:
            rkf,fkf = vals[inode]
            vs = ifoms(inode)
            rchi,rchib,rac = vs[0]
            fchi,fchib,fac = vs[1]
            #print ' inode chi ',inode,rchi,fchi
            #print ' inode chib ',inode,rchib,fchib
            #print ' inode ac   ',inode,rac,fac

            if (inode == 0): # inode == 0 is the vertex!!
                self.root.fill(self.prefix+'vtx_fchi',fchi)
                self.root.fill(self.prefix+'vtx_rchi',rchi)
                self.root.fill(self.prefix+'vtx_fchib',fchib)
                self.root.fill(self.prefix+'vtx_rchib',rchib)
                self.root.fill(self.prefix+'vtx_fac',fac)
                self.root.fill(self.prefix+'vtx_rac',rac)
                fnodes = fkf.nodes
                tchi2 = map(lambda nd: nd.getchi2(MTYPE),fnodes)
                for i in range(len(tchi2)): 
                    self.root.fill(self.prefix+'vtx_fchii_pf',i,tchi2[i])
                rnodes = rkf.nodes
                tchi2 = map(lambda nd: nd.getchi2(MTYPE),rnodes)
                for i in range(len(tchi2)): 
                    self.root.fill(self.prefix+'vtx_rchii_pf',i,tchi2[i])
            if (inode<=-5):
                self.root.fill(self.prefix+'nvtx_fchi',fchi)
                self.root.fill(self.prefix+'nvtx_fchib',fchib)
                self.root.fill(self.prefix+'nvtx_fac',fac)
                self.root.fill(self.prefix+'vnv_fchi',fchi,vfchi)
                self.root.fill(self.prefix+'vnv_fchib',fchib,vfchib)
                self.root.fill(self.prefix+'vnv_fac',fac,vfac)
                fnodes = fkf.nodes
                tchi2 = map(lambda nd: nd.getchi2(MTYPE),fnodes)
                for i in range(len(tchi2)): 
                    self.root.fill(self.prefix+'nvtx_fchii_pf',i,tchi2[i])
            if (inode>=5):
                self.root.fill(self.prefix+'nvtx_rchi',rchi)
                self.root.fill(self.prefix+'nvtx_rchib',rchib)
                self.root.fill(self.prefix+'nvtx_rac',rac)
                self.root.fill(self.prefix+'vnv_rchi',rchi,vrchi)
                self.root.fill(self.prefix+'vnv_rchib',rchib,vrchib)
                self.root.fill(self.prefix+'vnv_rac',rac,vrac)
                rnodes = rkf.nodes
                tchi2 = map(lambda nd: nd.getchi2(MTYPE),rnodes)
                for i in range(len(tchi2)): 
                    self.root.fill(self.prefix+'nvtx_rchii_pf',i,tchi2[i])

            self.root.fill(self.prefix+'faci',inode,fac)
            self.root.fill(self.prefix+'raci',inode,rac)
            self.root.fill(self.prefix+'aci',inode,fac+rac)
            self.root.fill(self.prefix+'fchii',inode,fchi)
            self.root.fill(self.prefix+'rchii',inode,rchi)
            self.root.fill(self.prefix+'chii',inode,0.5*(rchi+fchi))
            self.root.fill(self.prefix+'fchibi',inode,fchib)
            self.root.fill(self.prefix+'rchibi',inode,rchib)
            self.root.fill(self.prefix+'chibi',inode,0.5*(rchib+fchib))

            self.root.fill(self.prefix+'fac_pf',inode,fac)
            self.root.fill(self.prefix+'rac_pf',inode,rac)
            self.root.fill(self.prefix+'ac_pf',inode,fac+rac)
            self.root.fill(self.prefix+'fchi_pf',inode,fchi)
            self.root.fill(self.prefix+'rchi_pf',inode,rchi)
            self.root.fill(self.prefix+'chi_pf',inode,0.5*(rchi+fchi))
            self.root.fill(self.prefix+'fchib_pf',inode,fchib)
            self.root.fill(self.prefix+'rchib_pf',inode,rchib)
            self.root.fill(self.prefix+'chib_pf',inode,0.5*(rchib+fchib))

            def asym(a,b): 
                if (a+b) == 0.: return -1.
                return (a-b)/(a+b)

            self.root.fill(self.prefix+'v0a_fchii',inode,asym(fchi,vfchi))
            self.root.fill(self.prefix+'v0a_rchii',inode,asym(rchi,vrchi))
            self.root.fill(self.prefix+'v0a_fchibi',inode,asym(fchib,vfchib))
            self.root.fill(self.prefix+'v0a_rchibi',inode,asym(rchib,vrchib))

            self.root.fill(self.prefix+'v0a_fchi_pf',inode,asym(fchi,vfchi))
            self.root.fill(self.prefix+'v0a_rchi_pf',inode,asym(rchi,vrchi))
            self.root.fill(self.prefix+'v0a_fchib_pf',inode,asym(fchib,vfchib))
            self.root.fill(self.prefix+'v0a_rchib_pf',inode,asym(rchib,vrchib))

        return True


class HistosKFChiScan(IAlg):
    """ Alforithm to get the chi2-scan of the a kalman-filter
    """

    def define(self):
        self.prefix='chiscan_trk'
        self.nhits = 50
        self.ntracks = 2
        self.chimax = 30.
        self.path = 'rec/kfs'
        return

    def initialize(self):
        for i in range(self.ntracks):
            self.root.h2d(self.prefix+str(i)+'_fchii',self.nhits,0,self.nhits,40,0.,30.)
            self.root.h2d(self.prefix+str(i)+'_schii',self.nhits,0,self.nhits,40,0.,30.)
            self.root.hprf(self.prefix+str(i)+'_fchii_pf',self.nhits,0,self.nhits,0.,30.)
            self.root.hprf(self.prefix+str(i)+'_schii_pf',self.nhits,0,self.nhits,0.,30.)
        return True

    def execute(self):
        kfs = self.evt[self.path]
        if (not kfs): return False
        for i,kf in enumerate(kfs):
            fchis = map(lambda nd: nd.getchi2('filter'),kf.nodes)
            schis = map(lambda nd: nd.getchi2('smooth'),kf.nodes)
            for k,chi in enumerate(fchis):
                if (chi>self.chimax): continue
                self.root.fill(self.prefix+str(i)+'_fchii',k,chi)
                self.root.fill(self.prefix+str(i)+'_fchii_pf',k,chi)
            for k,chi in enumerate(schis):
                if (chi>self.chimax): continue
                self.root.fill(self.prefix+str(i)+'_schii',k,chi)
                self.root.fill(self.prefix+str(i)+'_schii_pf',k,chi)
        return True

#-------------------------------
# Event display
#-------------------------------

class EventDisplayPull(IAlg):
    """ Alforithm to event display the kalfman filter result
    """

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

class EventDisplayScan(IAlg):
    """ Algorithm to event display the scan of the track
    """

    def define(self):
        self.path = 'rec/kfs/bb/scan'
        self.prefix = 'evt_sn_'
        self.nnodes = 30
        self.nevents = 50
        self.ievt = 0 
        self.zsize = 12.
        return

    def execute(self):
        print " Scan! "
        vals = self.evt['rec/kfs/scan']
        if (not vals): return False
        inodes = vals.keys()
        inodes.sort()
        #print 'inodes ',inodes
        st=str(self.ievt)
        
        self.root.h1d(self.prefix+st+'_fchi',2*self.nnodes,-self.nnodes,self.nnodes)
        self.root.h1d(self.prefix+st+'_fchib',2*self.nnodes,-self.nnodes,self.nnodes)
        self.root.h1d(self.prefix+st+'_fac',2*self.nnodes,-self.nnodes,self.nnodes)
        self.root.h1d(self.prefix+st+'_rchi',2*self.nnodes,-self.nnodes,self.nnodes)
        self.root.h1d(self.prefix+st+'_rchib',2*self.nnodes,-self.nnodes,self.nnodes)
        self.root.h1d(self.prefix+st+'_rac',2*self.nnodes,-self.nnodes,self.nnodes)

        self.root.h2d(self.prefix+st+'_xyene',100,-self.zsize,self.zsize,100,-self.zsize,self.zsize)
        self.root.h2d(self.prefix+st+'_xyfchi',100,-self.zsize,self.zsize,100,-self.zsize,self.zsize)
        self.root.h2d(self.prefix+st+'_xyfchib',100,-self.zsize,self.zsize,100,-self.zsize,self.zsize)
        self.root.h2d(self.prefix+st+'_xyrchi',100,-self.zsize,self.zsize,100,-self.zsize,self.zsize)
        self.root.h2d(self.prefix+st+'_xyrchib',100,-self.zsize,self.zsize,100,-self.zsize,self.zsize)

        def ifoms(inode):
            rkf,fkf = vals[inode]
            rvals = 0.,0.,0.
            if (rkf): rvals = foms(rkf.nodes)
            fvals = 0.,0.,0.
            if (fkf): fvals = foms(fkf.nodes)
            return (rvals,fvals)

        vs = ifoms(0)
        vrchi,vrchib,vrac = vs[0]
        vfchi,vfchib,vfac = vs[1]

        for inode in inodes:
            rkf,fkf = vals[inode]
            if (not rkf or not fkf): continue
            vs = ifoms(inode)
            rchi,rchib,rac = vs[0]
            fchi,fchib,fac = vs[1]
            #print ' inode chi ',inode

            self.root.fill(self.prefix+st+'_fchi',inode,fchi)
            self.root.fill(self.prefix+st+'_fchib',inode,fchib)
            self.root.fill(self.prefix+st+'_fac',inode,fac)
            self.root.fill(self.prefix+st+'_rchi',inode,rchi)
            self.root.fill(self.prefix+st+'_rchib',inode,rchib)
            self.root.fill(self.prefix+st+'_rac',inode,rac)

            node0 = fkf.nodes[0]
            xx = node0.getstate('true').vec
            x,y,ene = xx[0],xx[1],xx[4]
            #print ' x, y, z, ene ',x,y,node0.zrun,ene
            self.root.fill(self.prefix+st+'_xyfchi',x,y,fchi)
            self.root.fill(self.prefix+st+'_xyfchib',x,y,fchib)
            self.root.fill(self.prefix+st+'_xyene',x,y,ene)
            self.root.fill(self.prefix+st+'_xyrchi',x,y,rchi)
            self.root.fill(self.prefix+st+'_xyrchib',x,y,rchib)

        names = map(lambda x: self.prefix+st+'_'+x,['fchi','fchib','fac','rchi','rchib','rac'])
        hs = map(lambda name: self.root.get(name),names)
        t0 = tcanvas(hs,nx=3,name='evtscan_'+st)
        raw_input('press key ')    

        names = map(lambda x: self.prefix+st+'_'+x,['xyfchi','xyfchib','xyene','xyrchi','xyrchib'])
        hs = map(lambda name: self.root.get(name),names)
        t0 = tcanvas(hs,nx=2,name='evtscan_'+st)
        raw_input('press key ')    

        self.ievt+=1

        return True

class EventDisplay(IAlg):
    """ Algorithm to event display of the generated events
    """

    def define(self):
        self.path = 'sim/zstates'
        self.prefix = 'evt_'
        self.nevents = 50
        self.zsize = 15.
        self.ievt = 0 
        return

    def execute(self):
        #print " Scan! "
        zsegs = self.evt[self.path]
        if (not zsegs): return False
        st=str(self.ievt)
        #print " Scan! ",st

        self.root.h2d(self.prefix+st+'_xy',100,-self.zsize,self.zsize,100,-self.zsize,self.zsize)
        self.root.h2d(self.prefix+st+'_xz',100,-self.zsize,self.zsize,100,-self.zsize,self.zsize)
        self.root.h2d(self.prefix+st+'_yz',100,-self.zsize,self.zsize,100,-self.zsize,self.zsize)

        for zseg in zsegs:
            for zst in zseg:    
                x,y,z,ux,uy,uz,ee = zst
                self.root.fill(self.prefix+st+'_xy',x,y,ee)
                self.root.fill(self.prefix+st+'_xz',x,z,ee)
                self.root.fill(self.prefix+st+'_yz',y,z,ee)
                
        names = map(lambda x: self.prefix+st+'_'+x,['xy','xz','yz'])
        hs = map(lambda name: self.root.get(name),names)

        t0 = tcanvas(hs,nx=2,name='evt_'+st)
        raw_input('press key ')    
        self.ievt+=1

        return True

#-----------------------------
# Alex programs
#-----------------------------

def doscan():
    """ Alex program to generate doble beta events and do a scan along the track
    """

    # create application
    alex = Alex()
    alex.nevts = 2000
    root = ROOTSvc('root','akfbeta.root')
    alex.addsvc(root)

    # simulation
    #----------------

    # generate double beta
    agenstates = GenerateDoubleBeta('agenstates')
    #agenstates = GenerateBeta('agenstates')
    alex.addalg(agenstates)

    # histogram simulated states
    hisgenstates = HistosGenerateStates('hisgenstates')  
    hisgenstates.imports.append('root')
    alex.addalg(hisgenstates)

    # genrate digits
    agendigits = CreateDigits('agendigits')
    agendigits.dz = 0.5
    alex.addalg(agendigits)  

    # provide a track with the two electrons
    betanodes = DoubleBetaNodes('betanodes')
    alex.addalg(betanodes)

    # fitting
    #------------

    # do forward fit
    akf = FitNodes('kfbetanodes')
    akf.path = 'rec/nodes/bb'
    akf.opath = 'rec/kfs/bb'
    alex.addalg(akf)

    # histo forward fit
    hisakf = HistosKFFit('hisakf')
    hisakf.imports.append('root')
    hisakf.prefix='kf_'
    hisakf.path='rec/kfs/bb'
    alex.addalg(hisakf)

    # histo chi2-scan of the forward fit
    hischiscan = HistosKFChiScan('hischiscan')
    hischiscan.imports.append('root')
    hischiscan.path = 'rec/kfs/bb'
    alex.addalg(hischiscan)

    # reverse the nodes to do the reverse fit
    arevnodes = ReverseNodes('arevnodes')
    arevnodes.path = 'rec/nodes/bb'
    arevnodes.opath = 'rec/nodes/bb/rev'
    alex.addalg(arevnodes)

    # do the reverse fit
    akfrev = FitNodes('akfrev')
    akfrev.path = 'rec/nodes/bb/rev'
    akfrev.opath = 'rec/kfs/bb/rev'
    alex.addalg(akfrev)

    # histotgram the reverse fit
    hisakfrev = HistosKFFit('hisakfrev')
    hisakfrev.imports.append('root')
    hisakfrev.prefix='kfrev_'
    hisakfrev.path='rec/kfs/bb/rev'
    alex.addalg(hisakfrev)

    # histogram the chi2-scan of the reverse fit
    hischiscanrev = HistosKFChiScan('hischiscanrev')
    hischiscanrev.imports.append('root')
    hischiscanrev.path='rec/kfs/bb/rev'
    hischiscanrev.prefix='chiscanrev_trk'
    alex.addalg(hischiscanrev)

    # scan along the track
    #----------------------
   
    # do the analysis of the forward/reverse fit
    anarev = HistosAnaReverse('anrev')
    anarev.imports.append('root')
    anarev.pathnor = 'rec/kfs/bb'
    anarev.pathrev = 'rec/kfs/bb/rev'
    alex.addalg(anarev)

    akf = FitScanNodes('ascan')
    akf.path = 'rec/nodes/bb'
    akf.opath = 'rec/kfs/bb/scan'
    akf.iscan = 4
    alex.addalg(akf)

    hakf = HistosScanNodes('hscan')
    hakf.imports.append('root')
    alex.addalg(hakf)

    # some event diplays
    #---------------

    #evtakf = EventDisplayScan('evtakf')
    #evtakf.imports.append('root')
    #alex.addalg(evtakf)

    #evtdisplay = EventDisplay('evtdisplay')
    #evtdisplay.imports.append('root')
    #alex.addalg(evtdisplay)

    alex.run()
    
def dotracks():
    """ Alex program to generate a single electron or doble beta event with toy-MC.
    and fit in forward and reverse mode. Histograms are produced.
    """

    # create aplication
    alex = Alex()
    alex.msg.level = Message.Info
    alex.nevts = 1
    root = ROOTSvc('root','akfbeta.root')
    alex.addsvc(root)

    # simulate states
    #--------------------------
    agenstates = GenerateDoubleBeta('agenstates') # uncoment to generate doble-beta
    #agenstates = GenerateBeta('agenstates')
    agenstates.E0 = 2.5
    alex.addalg(agenstates)

    # histogram simulated states
    hisgenstates = HistosGenerateStates('hisgenstates')  
    hisgenstates.imports.append('root')
    alex.addalg(hisgenstates)

    # generate hits
    agendigits = CreateDigits('agendigits')
    agendigits.dz = 0.5
    alex.addalg(agendigits)  

    # create nodes
    #agennodes = CreateNodes('agennodes')
    #alex.addalg(agennodes)

    # uncoment to generetate doble-beta nodes in one track
    betanodes = DoubleBetaNodes('betanodes')
    alex.addalg(betanodes)

    # fit
    #-----------------

    # fit in the forward direction
    akf = FitNodes('akf')
    akf.path = 'rec/nodes'
    akf.opath = 'rec/kfs'
    alex.addalg(akf)

    # histogram the forward fit
    hisakf = HistosKFFit('hisakf')
    hisakf.imports.append('root')
    hisakf.prefix='kf_'
    hisakf.path='rec/kfs'
    alex.addalg(hisakf)

    # histogram the chi2-scan of the forward fit
    hischiscan = HistosKFChiScan('hischiscan')
    hischiscan.imports.append('root')
    hischiscan.path = 'rec/kfs'
    alex.addalg(hischiscan)

    # reverse the nodes to do a reverse fit
    arevnodes = ReverseNodes('arevnodes')
    arevnodes.path = 'rec/nodes'
    arevnodes.opath = 'rec/nodes/rev'
    alex.addalg(arevnodes)

    # do reverse fit
    akfrev = FitNodes('akfrev')
    akfrev.path = 'rec/nodes/rev'
    akfrev.opath = 'rec/kfs/rev'
    alex.addalg(akfrev)

    # histogram the reverse fit
    hisakfrev = HistosKFFit('hisakfrev')
    hisakfrev.imports.append('root')
    hisakfrev.prefix='kfrev_'
    hisakfrev.path='rec/kfs/rev'
    alex.addalg(hisakfrev)

    # histogram the chi2-scan of the reverse fit
    hischiscanrev = HistosKFChiScan('hischiscanrev')
    hischiscanrev.imports.append('root')
    hischiscanrev.path='rec/kfs/rev'
    hischiscanrev.prefix='chiscanrev_trk'
    alex.addalg(hischiscanrev)

    # Compare forward/reverse fit
    #-----------------------------

    # do the analysis of the forward/reverse fit
    anarev = HistosAnaReverse('anrev')
    anarev.imports.append('root')
    alex.addalg(anarev)

    alex.run()    
    
if __name__ == '__main__':
    """ Alex program: two configurations:
    doscan - scan a track to find the vertex
    dotrack - fit the track in forward/reverse mode
    uncomment the configuration to run.
    """
    
    #doscan()
    dotracks()
