"""
trackgen.py

Toy MC generates tracks with multiple scattering.
Tracks are stores in a ROOT tree

Tracks can be recover from the Tree.

Track is a list of states
Each State has
    x y z ux uy uz E

    x,y,z point of the track
    ux,uy,uz: direction vector 
    E: energy of the track 
    
    All distances are in cm, Energy in MeV?
"""

from math import *
from alex.ialex import IAlg
from alex.alex import Alex
from alex.reader import TreeReader, Reader
from alex.rootsvc import ROOTSvc, ttree, h1d, h2d
from alex.algorithms import Dump
import agenparam
import random
import scipy.integrate as integrate
from ROOT import gDirectory

#-----------------
#   ToyMC code
#-----------------

def momentum(ene,mass = agenparam.masse):
    """ return the momentum from a given energy 
    """
    p = sqrt((ene+mass)**2-mass**2)
    return p

def theta_ms(p,dis,mass=agenparam.masse):
    """ return the theta MS angle for a particle with mass, momentum (p) that travel a distance (dis)
    """
    ene = sqrt(mass**2+p**2)
    beta = p/ene
    tms = (13.6)/(p*1.*beta)*sqrt(dis*1.)*(1+0.038*log(dis*1.))
    return tms

class Track:
    """ class Tracks, holds a list of states
    Each state has (x,y,z,ux,uy,uz,ee), the position, the steer-cos and the energy
    """

    def __init__(self,mass=agenparam.masse,states=[]):
        self.mass = mass
        self.states = list(states)
        return        

    def __str__(self):
        return str(self.states)

    def __repr__(self):
        return str(self.states)

def genstate(state):
    """ generate a next state from this state (state), the particile losses deltaE in the step
    """
    xi,yi,zi,uxi,uyi,uzi,ei = state
    #self.msg.verbose("# state -ini ",state)
    pi = momentum(ei)
    deltaE = min(agenparam.eslice,ei)
    ef = max(ei-deltaE,0.)
    deltaX = integrate.quad(lambda ee: 1./agenparam.xesp(ee),
                            ef,ei,limit=1000)[0]
    dx = deltaX*uxi; dy = deltaX*uyi; dz = deltaX*uzi;
    xf = xi+dx; yf = yi+dy; zf = zi+dz;
    #self.msg.verbose('#  x,y,z,delta ',xf,yf,zf,deltaX)

    # Determine the scattering angles in the frame in which the track
    #  direction is aligned with the z-axis.
    x0 = deltaX/agenparam.Lr
    x0 = min(0.001,max(100.,x0))
    sigma_theta = theta_ms(pi,x0);
    
    thetaX = random.gauss(0,sigma_theta);
    thetaY = random.gauss(0,sigma_theta);
    tanX = tan(thetaX);
    tanY = tan(thetaY);
    #self.msg.verbose('#  ms theta, tanX, tanY ',sigma_theta,tanX,tanY)
    
    # Compute the direction cosines of the rotation matrix to move to the lab frame.
    nxy = sqrt(uxi**2 + uyi**2);
    if(nxy > 0.):
        alpha1 = uyi/nxy; alpha2 = -uxi*uzi/nxy; alpha3 = uxi;
        beta1 = -uxi/nxy; beta2 = -uyi*uzi/nxy; beta3 = uyi;
        gamma1 = 0.; gamma2 = nxy; gamma3 = uzi;
    else:
        # Special case; the direction vector is the x-axis.  Choose
        #  the orthonormal basis as the normal unit vectors.
        alpha1 = 1.; alpha2 = 0.; alpha3 = 0.;
        beta1 = 0.; beta2 = 1.; beta3 = 0.;
        gamma1 = 0.; gamma2 = 0.; gamma3 = 1.;

    # Determine direction vector components in the reference (lab) frame.
    nrm = sqrt(tanX**2 + tanY**2 + 1);
    xp = (alpha1*tanX + alpha2*tanY + alpha3)/nrm;
    yp = (beta1*tanX + beta2*tanY + beta3)/nrm;
    zp = (gamma1*tanX + gamma2*tanY + gamma3)/nrm;
    
    # Set the new direction vector.
    nrm = sqrt(xp**2 + yp**2 + zp**2);
    uxf = xp/nrm;
    uyf = yp/nrm;
    uzf = zp/nrm;

    #self.msg.verbose('#  next ux,uy,uz ',uxf,uyf,uzf)
    
    state = (xf,yf,zf,uxf,uyf,uzf,ef)
    #self.msg.verbose('# state -next ',state)
    return state

def gentrack():
    """ generate a electron track with initial (0.,0.,0.,0.,0.,1.,E0)
    """
    track = Track()
    state0 = (0.,0.,0.,0.,0.,1.,agenparam.E_0)
    ene = state0[-1]
    while (ene>0.):
        state = genstate(state0)
        track.states.append(state)
        state0=state
        ene = state0[-1]
    return track

#--------------------
#     Alex code
#--------------------

class TrackGen(IAlg):
    """ IAlg to generate a track and store it in a ROOT Tree
    """
 
    def define(self):
        """ path : path to store the track in evtstore
        """
        self.path = self.name 
        return

    def initialize(self):
        """ book the ntuple and put it in ROOTSvc
        """
        labels = [('x','F',51),('y','F',51),('z','F',51),
                  ('ux','F',51),('uy','F',51),('uz','F',51),
                  ('ee','F',51)]
        tree = ttree(self.name,labels)
        self.root.put(tree)
        return True

    def execute(self):
        """ create the track, store it in the evtstore and fill the ntuple
        """
        track = gentrack()
        self.evt[self.name]=track
        x = map(lambda s: s[0],track.states)
        y = map(lambda s: s[1],track.states)
        z = map(lambda s: s[2],track.states)
        ux = map(lambda s: s[3],track.states)
        uy = map(lambda s: s[4],track.states)
        uz = map(lambda s: s[5],track.states)
        ee = map(lambda s: s[6],track.states)
        tup = {'x':x,'y':y,'z':z,'ux':ux,'uy':uy,'uz':uz,'ee':ee}
        self.root.fill(self.name,tup)
        return True        

class TrackRecover(IAlg):
    """ IAlg to recover the track from the read data and store it in the evtstore
    """
    
    def define(self):
        """ ipath: path from where to recover data
            opath: path to store the rocovered track
        """
        self.ipath = 'reader' #input path of the data
        self.opath = self.name #output path for the recover track
        return

    def execute(self):
        """ recover the track from the data in the store
        
        """
        data = self.evt.get(self.ipath)
        if (not data): return False
        n = len(data['ee'])
        track = Track()
        for i in range(n):
            state = (data['x'][i],data['y'][i],data['z'][i],
                     data['ux'][i],data['uy'][i],data['uz'][i],data['ee'][i])
            self.msg.verbose(' state ',state)
            track.states.append(state)
        self.evt[self.opath] = track
        return True

class TrackHistos(IAlg):

    def define(self):
        self.path = 'track'
        return
    
    def initialize(self):
        self.root.h1d('hx',100,-30.,30.)
        self.root.h1d('hy',100,-30.,30.)
        self.root.h1d('hz',100,0.,40.)
        self.root.h1d('hux',100,-1.2,1.2)
        self.root.h1d('huy',100,-1.2,1.2)
        self.root.h1d('huz',100,-1.2,1.2)
        self.root.h1d('hee',100,0.,3.)
        return True

    def execute(self):
        track = self.evt.get(self.path)
        if (not track): return False
        for state in track.states:
            x,y,z,ux,uy,uz,ee = state
            self.root.fill('hx',x)
            self.root.fill('hy',y)
            self.root.fill('hz',z)
            self.root.fill('hux',ux)
            self.root.fill('huy',ux)
            self.root.fill('huz',uz)
            self.root.fill('hee',ee)
            self.root.fill('hxz',z,x)
            self.root.fill('hyz',z,y)
            self.root.fill('huxz',z,ux)
            self.root.fill('huyz',z,uy)
            self.root.fill('hez',z,ee)
        return True
            

class HistoDir(IAlg):

    def execute(self):
        keys = self.root.tstore.keys()
        for key in keys:
            obj = self.root.fetch(key)
            print '# histo ',key,' dir ',obj.GetDirectory(),' file',obj.GetDirectory().GetFile()
        return True
#---------------------------
#    Alex Jobs
#---------------------------


def jobtrackgen():
    """ Job to generate ntracks and store them into a ntuple
    """
    
    # confifure aplication 
    alex = Alex('alex')
    alex.nevts = 100
    
    root = ROOTSvc('root',fname='atrackgen.root')
    alex.addsvc(root)
    
    # configure flow
    gen = TrackGen('trackgen')
    gen.imports.append('root')
    alex.addalg(gen)
    
    histos = TrackHistos('histos')
    histos.imports.append('root')
    histos.path = 'trackgen'
    alex.addalg(histos)

    root.h2d('hxz',100,0.,40.,100,-30.,30.)
    root.h2d('hyz',100,0.,40.,100,-30.,30.)
    root.h2d('huxz',100,0.,40.,100,-1.2,1.2)
    root.h2d('huyz',100,0.,40.,100,-1.2,1.2)
    root.h2d('hez',100,0.,40.,100,0.,2.6)

    # run!
    alex.run()
    return
    
def jobtrackrecover():
    """ Job to recover ntracks and fill some histograms
    """    
    
    # confifure aplication 
    alex = Alex('alex')
    alex.nevts = 400
    
    root = ROOTSvc(name='root',fname='atrackrecover.root')
    alex.addsvc(root)

    # configure flow
    RR = lambda fname : TreeReader(fname,tname='trackgen')
    reader = Reader(name='reader',freader=RR)
    #read.fnames = ['trackgen1.root','trackgen2.root','trackgen3.root']
    reader.fnames = ['atrackgen.root']
    alex.addalg(reader)

    recover = TrackRecover('recover')
    alex.addalg(recover)
        
    histos = TrackHistos('histos')
    histos.imports.append('root')
    histos.path = 'recover'
    alex.addalg(histos)

    root.h2d('hxz',100,0.,40.,100,-30.,30.)
    root.h2d('hyz',100,0.,40.,100,-30.,30.)
    root.h2d('huxz',100,0.,40.,100,-1.1,1.1)
    root.h2d('huyz',100,0.,40.,100,-1.1,1.1)
    root.h2d('hez',100,0.,40.,100,0.,3.)

    # run!
    alex.run()
    return


if __name__ == '__main__':
    #jobtrackgen()
    jobtrackrecover()
