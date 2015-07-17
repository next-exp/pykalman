"""

Physics data and formulae for NEXT

"""
from KFBase import KFVector
from KFBase import KFMatrix, KFMatrixNull, KFMatrixUnitary
from KFBase import Random
from math import *
import random
import numpy as np
from scipy.interpolate import interp1d as sc_interpol
from scipy.optimize import bisect as sc_root
from scipy.integrate import quad as sc_integrate
from exceptions import ZeroDivisionError

DEBUG = False
WARNING = True

MAXTHETA = float('inf')#1e6*pi
ENEMIN = 0.0

def debug(comment,arg=''):
    if (DEBUG): print comment,arg

def warning(comment,arg=''):
    if (WARNING): print comment,arg

def energy(p,m=0.511):
    ene = sqrt(p*p+m*m)
    return ene

def kinenergy(p,m=0.511):
    """ return the kinetic energy for a momentum p """
    ene = energy(p,m)
    te = ene-m
    return te

# Be careful with energy, kin-energy and p
# MS works with p and energy
# Eloss works with kin-energy
# Qbb works with kin-energy    

def thetams(p,dis,X0,mass=0.511):
    """ return the theta MS angle for a particle with mass, momentum (p) 
    that travel a distance (dis) in a radiation length X0
    """
    dis = abs(dis)
    if (dis == 0.): return 0.
    if (X0 <=0.): return 0.
    udis =dis/X0
    ene = sqrt(mass**2+p**2)
    beta = p/ene
    tms = (13.6)/(p*1.*beta)*sqrt(udis*1.)*(1+0.038*log(udis*1.))
    debug('thetams p, udis, theta ',(p,udis,tms))
    return tms

def UMatrix(udir):
    """ returns the rotation matrix that transfrom vectors in the track ref. system to the global system
    """
    ux,uy,uz = udir
    uu = sqrt(ux*ux+uy*uy+uz*uz)
    ux = ux/uu; uy = uy/uu; uz = uz/uu
    #print ' u ',ux,uy,uz
    if (abs(uz)==0.): raise ZeroDivisonError
    vx = uz; vy = 0.; vz = -ux
    vv = sqrt(vx*vx+vy*vy+vz*vz)
    if (abs(vv)==0.): raise ZeroDivisonError
    vx = vx/vv; vy = vy/vv; vz=vz/vv
    #print ' v ',vx,vy,vz
    wx = -uy*ux/vv; wy = vv; wz = -uy*uz/vv
    ww = sqrt(wx*wx+wy*wy+wz*wz)
    wx = wx/ww; wy = wy/ww; wz = wz/ww        
    #print ' w ',wx,wy,wz
    # UMatrix convert track ref system to global system
    U = KFMatrix( [[vx,wx,ux ],[vy,wy,uy],[vz,wz,uz]])
    debug('UMatrix udir, U ',(udir,U))
    return U    

class MSNoise:
    """ Multiple scattering helper class
    it provides thetams, cov, random values 
    it takes as input a radiation length X0, and a mass
    """
    
    def __init__(self,X0,mass=0.511,thetamax=MAXTHETA):
        """ constructor of a multiple scattering noiser with a X0
        """
        self.X0 = X0
        self.mass = mass
        self.thetamax = thetamax
        return

    def validstep(self,p,dis):
        """ return true is this step is valid (if thetams is not so large!)
        """
        tms = self.theta(p,dis)
        ok = (tms<=self.thetamax)
        if (not ok):
            warning("msnoiser.validstep p,dis,tms,ok ",(p,dis,tms,ok))
        debug("msnoiser.validstep p,dis,tms,ok ",(p,dis,tms,ok))
        return ok

    def theta(self,p,dis):
        """ return theta of ms
        """
        tms = thetams(p,dis,self.X0,self.mass)
        debug('msnoiser.theta p, dis  tms',tms)
        return tms
 
    def Q0Matrix(self,p,dis):
        """ returns the multiple scattering Cov matrix in the track system
        (z-direction if the direction of the track).
        particle has momentum p and travels a distance c
        (1,2) coordinates are transverse coordinates respect track direction
        """
        theta0 = self.theta(p,dis)
        theta02 = theta0*theta0
        Q = KFMatrixNull(4,4)
        Q[0,0] = dis*dis*theta02/3.
        Q[0,2] = dis*theta02/2.
        Q[1,1] = dis*dis*theta02/3.
        Q[1,3] = dis*theta02/2.
        Q[2,2] = theta02
        Q[3,3] = theta02
        for i in range(4):
            for j in range(i+1,4):
                Q[j,i] = Q[i,j]
        debug('msnoiser.Q0Matrix p,dis,Q ',(p,dis,Q))
        return Q

    def random(self,p,dis):
        """ generate multiple scattering random variables in x,theta.
        x is the transverse direction
        and theta the angle.
        """
        theta0 = MS.theta(p,dis)
        rho = sqrt(3.)/2.
        z1,z2 = random.gauss(0.,1.),random.gauss(0.,1.)
        y = z1*dis*theta0/sqrt(12.)+z2*dis*theta0/2.
        theta = z2*theta0
        debug('msnoise.random p,dis,x,theta ',(p,dis,y,theta))
        return y,theta
    
    def XUrandom(self,p,dis,x0,udir):
        """ return a position and direction (in the global system) 
        after a random MS of a particle with momentum p
        that traverses a distance dis with a direction udir
        """
        udir.Unit()
        U = UMatrix(udir)
        #print ' U ',U
        Q0 = self.Q0Matrix(p,dis)
        #print ' Q0 ',Q0
        x1,x2,t1,t2 = Random.cov(Q0)
        #print ' x1,x2,t1,t2 ',x1,x2,t1,t2
        t1 = tan(t1); t2 = tan(t2); nor = sqrt(1.+t1*t1+t2*t2)
        xt = KFVector([x1,x2,0.]) 
        ut = KFVector([t1/nor,t2/nor,1./nor]) 
        #print ' xt ',xt
        #print ' ut ',ut
        xf = x0 + dis*udir
        #print ' xf ',xf
        xf = xf+U*xt
        #print ' xf ',xf
        uf = U*ut
        #print ' uf ',uf
        if (DEBUG): 
            print 'msnoise.XUrandom p,dis,x0,udir ',p,dis,x0,udir
            print 'msnoise.XUrandom xt,ut ',xt,ut
            print 'msnoise.XUrandom xf,uf ',xf,uf
        return xf,uf

class ELoss:
    """ empty class for energy loss
    """
    
    mindeltax = 0.2
    
    def validstep(self,ene,ds):
        """ return true is this is a valid step for this energy and distance
        """
        return True
        
    def deltae(self,ene0,deltax=0.):
        """ return delta-e associated to a delta-x distance
        """
        # temporal
        return 0.,deltax

    def deltax(self,ene0,deltaene=0.):
        """ return delta-x associated to a deltaene loss
        """
        return 0.,1.

class HPXeELoss(ELoss):
    """ Energy loss for electrons in High Pressure Xenon
    """

    def __init__(self,rho,enemin=ENEMIN):
        """ construction of HPXe Energy loss for election.
        it reuires the rho, Xe density
        """
        self.rho = rho
        self.enemin=enemin
        # Read in the stopping power and interpolate.
        xesp_tbl = np.loadtxt("/home/brais/Documents/Next/pykalman/toyMC/xe_estopping_power_NIST.dat");
        self.evals = xesp_tbl[:,0]
        self.dEdxvals = xesp_tbl[:,1];
        self.evals = np.insert(self.evals,0,0.0);
        self.dEdxvals = np.insert(self.dEdxvals,0,self.dEdxvals[0]);
        self.xesp = sc_interpol(self.evals,self.dEdxvals*self.rho,kind='cubic')
        # this is the dE/dx(E) function
        return 

    def validstep(self,b,ds):
        """ return true is this is a valid step for this energy and distance
        """
        ok = True
        ds = abs(ds)
        de,dz = self.deltae(ene,ds)
        if (dz+0.001<ds): 
            ok = False
            print '1' #   ;
        if (ene<self.enemin): 
            ok = False
            print '2'#    ;
        if (ene-de<=0.): 
            ok= False
            print '3' #;
        if (not ok):
            debug("HPXeELoss.validstep ene0,de,ds,ok ",(ene,de,ds,ok))
        return ok
        
    def deltae(self,ene0,deltax=0.01):
        """ return delta-e associated to a delta-x distance
        """
        # temporal @TODO revisit this
        deltax = abs(deltax)
        de = 0.03556*(deltax/0.5)
        de = min(ene0,de)
        deltax = (de/0.03556)*0.5
        #ff = lambda dene : deltax - self.deltax(ene0,dene)[1]
        #if (ff(0.)*ff(ene0)>0.): 
        #    return ene0,self.deltax(ene0,ene0)
        #de = sc_root(ff,0.,ene0) # find the root
        debug('HPXeELoss.deltae ene,de,ds ',(ene0,de,deltax))
        return de,deltax

    def deltax(self,ene0,deltaene=0.01):
        """ return delta-x associated to a deltaene loss
        """
        # Determine the distance of this step.
        deltaene=abs(deltaene)
        enef = max(0.,ene0-deltaene)
        de = ene0-enef
        #dis = 0.5*(de/0.03556)
        dis = sc_integrate(lambda ee: 1./self.xesp(ee),enef,ene0,limit=1000)[0]
        debug('HPXeELoss.deltax ene,de,ds ',(ene0,de,dis))
        return de,dis
