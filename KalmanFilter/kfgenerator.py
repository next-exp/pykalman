from KFBase import KFVector 
from KFBase import KFMatrix,KFMatrixNull,KFMatrixUnitary,Random 
from kffilter import KFData,KFNode,KFModel,KFFilter
from math import *
import random

"""
KalmanFilter implementation for NEXT

The model is a stright line with (x,y,tx,ty,p) states

It has MS noise and energy loss

"""

DEBUG = False
WARNING = True

def debug(comment,arg=''):
    if (DEBUG): print "INFO ",comment,arg

def warning(comment,arg=''):
    if (WARNING): print "WARNING ",comment,arg

def energy(p,m=0.511):
    ene = sqrt(p*p+m*m)
    return ene

def kinenergy(p,m=0.511):
    """ return the kinetic energy for a momentum p """
    ene = energy(p,m)
    te = ene-m
    return te

def kinmomentum(kene,m=0.511):
    ene = kene+m
    p = sqrt(ene*ene-m*m)
    return p

class KFGenerator:
    """ Generator of particles states, a (x,y,z,ux,uy,ux,ene) list, 
    starts from a particle state seed and makes steps of energy loss
    it uses an msnoiser and eloss classes
    """

    def __init__(self,msnoiser,eloss,deltae=0.02,emin=0.05):
        """ constructor of the generator with a noiser, eloss and a range for delta energy ad distance
        """
        self.msnoiser = msnoiser
        self.eloss = eloss
        self.deltae = deltae
        self.emin = emin 
        return

    def step(self,state0):
        x0,y0,z0,ux0,uy0,uz0,ee0 = state0
        de,ds = self.eloss.deltax(ee0,self.deltae)
        ok = self.msnoiser.validstep(ee0,ds)
        if (not ok): return ok,state0
        xv0 = KFVector([x0,y0,z0])
        uv0 = KFVector([ux0,uy0,uz0])
        uv0.Unit()
        #print ' x0 ',xv0
        #print ' u0 ',uv0
        #print ' ee0 ',ee0
        #print ' de, ds ',de,ds
        pp = kinmomentum(ee0)
        xvf,uvf = self.msnoiser.XUrandom(pp,ds,xv0,uv0)
        x,y,z=xvf; ux,uy,uz=uvf
        ee = ee0-de
        state = [x,y,z,ux,uy,uz,ee]
        debug('kgenerator.step state ',state)
        return ok,state

    def generate(self,state0):
        """ generate nstates starting from state0, each one with an delta e loss
        """
        state = list(state0)
        states = [state]
        ok,ee = True,state[-1]
        while (ok and ee>self.emin):
            ok,state = self.step(state0)
            if (ok): states.append(state)
            state0 = list(state)
            ee = state[-1]
        debug('kfgenerator.generate states ',len(states))
        return states

def stateatz(states,z):
    def getzi(st0,st1):
        x0,y0,z0,ux0,uy0,uz0,ee0=st0
        x1,y1,z1,ux1,uy1,uz1,ee1=st1
        ok = False
        if (z0<=z and z<z1): ok =True
        if (z0>=z and z>z1): ok = True
        if (ok): return z
        return None
    def zstate(st0,st1,z):
        x0,y0,z0,ux0,uy0,uz0,ee0=st0
        x1,y1,z1,ux1,uy1,uz1,ee1=st1
        if (z-z0 == 0.): return st0
        dz = (z-z0)/(z1-z0)
        x = x0+(x1-x0)*dz
        y = y0+(y1-y0)*dz
        ux = ux0+(ux1-ux0)*dz
        uy = uy0+(uy1-uy0)*dz
        uz = uz0+(uz1-uz0)*dz
        uu = sqrt(ux*ux+uy*uy+uz*uz)
        ux=ux/uu;uy=uy/uu;uz=uz/uu
        ee = ee0-(ee0-ee1)*(z-z0)/(z1-z0)
        dene = ee-ee0
        stz = [x,y,z,ux,uy,uz,ee]
        #print " st0 ",st0
        #print " st1 ",st1
        #print " stz ",stz
        return stz
    nn = len(states)
    for ii in range(nn-1):
        st0,st1 = states[ii],states[ii+1]
        zi = getzi(st0,st1)
        if (not zi): continue
        zst =  zstate(st0,st1,zi)
    debug("kfgenerator.stateatz state ",zst)
    return zst


def zrunsample(states,zs,epsilon=0.01):
    """ sample the states, return states at zs positions
    """
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
        ux = ux0+ux1*(z-z0)/(z1-z0)
        uy = uy0+uy1*(z-z0)/(z1-z0)
        uz = uz0+uz1*(z-z0)/(z1-z0)
        uu = sqrt(ux*ux+uy*uy+uz*uz)
        ux=ux/uu;uy=uy/uu;uz=uz/uu
        ee = ee0-(ee0-ee1)*(z-z0)/(z1-z0)
        dene = ee-ee0
        stz = (x,y,z,ux,uy,uz,ee)
        #print " st0 ",st0
        #print " st1 ",st1
        #print " stz ",stz
        return stz
    udir = +1.
    for ii in range(0,nn-1):
        st0,st1 = states[ii],states[ii+1]
        zi = getzi(st0,st1)
        if (not zi): continue
        zst =  zstate(st0,st1,zi)
        zst = list(zst)
        if (len(zstates)>0): 
            zp = zstates[-1]
            zzp,zzn = zp[2],zst[2]
            if (zzn>zzp): udir=+1.
            elif (zzn<zzp): udir = -1.
            elif (zzn==zzp): zst[2]=zst[2]+udir*epsilon
        zstates.append(zst)
    debug("kfgenerator.zsample zstates ",zstates)
    return zstates
    
def zransample(states,p=0.2):
    """ sample random the states with a given probability
    """
    zs = []
    for state in states:
        pi = random.uniform(0.,1.)
        if (pi<=p): zs.append(state)
    debug("kfgenerator.zrunsample zstates ",zstates)
    return zs

def zavesample_old(states,n=4):
    """ sample the states making the average
    """
    zs = []
    m = int(len(states)/n)
    for i in range(m):
        st0 = [0,0,0,0,0,0,0]
        for k in range(n):
            kk = n*i+k
            for ii in range(7):
                st0[ii] += states[kk][ii]
        st0 = map(lambda xi: xi/n,st0)
        zs.append(st0)
    debug("kfgenerator.zavesample zstates ",zs)
    return zs
    
def zavesample(states,n=4):
    z0,zf = states[0],states[-1]
    states = states[1:-1]
    zsts = [z0]
    ns = len(states)
    mm = int(ns/n)
    for i in range(mm):
        istates = states[i*n:min((i+1)*n,ns)]
        iz = sum(map(lambda st: st[2],istates))/len(istates)
        izst = stateatz(istates,iz)
        #print 'zavesamples2 i states \n'
        #for state in istates: print '\t ',state
        #print 'zavesamples2 i zi ',iz
        #print 'zavesmaples2 i states ',izst
        zsts.append(izst)
    zsts.append(zf)
    debug('kfgenerator.zavesamples zsts ',zsts)
    #print 'zavesmaples2 states ',zsts
    return zsts

def zsegments(states):
    """ segment the states in forward-backward segments
    """
    segs = []
    def iturn(states):
        udir = states[0][5]
        for i in range(len(states)):
            fx = states[i][5]*udir
            if (fx<0): return i
        return i
    while (len(states)>2):
        i = iturn(states)
        seg = states[:i]
        segs.append(seg)
        states = states[i:]
    debug('kfgenerator.zsegments ',(len(segs),segs))
    return segs

def zrunsegments(nodes):
    """ segment the hits, states, nodes according with z
    """
    segs = []  
    def iturn(nodes):
        z0 = nodes[0].zrun
        zi = nodes[1].zrun
        udir = (zi-z0)
        if (udir == 0): return 1
        udir = udir/abs(zi-z0)
        for i in range(1,len(nodes)):
            zi = nodes[i].zrun
            d = (zi-z0)*udir
            if (d<=0.): return i
            z0 = zi
        return i
    while (len(nodes)>2):
        i = iturn(nodes)
        seg = nodes[:i]
        segs.append(seg)
        nodes = nodes[i:]
    #for seg in segs:
    #    print ' zrunsegments ',map(lambda nd: nd.zrun,seg)
    debug('kfzline.zrunsegments ',(len(segs),segs))
    return segs   

def zedeltas(enes):
    denes = []
    for i in enumerate(enes[:-1]):
        dene = enes[i+1]-ene[i]
    denes.append(enes[-1])
    return denes
