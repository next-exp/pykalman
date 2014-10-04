"""

Physics data and formulae for NEXT

"""

from KFBase import KFVector
from KFBase import KFMatrix, KFMatrixNull, KFMatrixUnitary
from math import *
import random
import numpy as np
from scipy.interpolate import interp1d

class Phys:

    # density of ideal gas at T=0C, P=1 atm in cm^(-3)
    gas_rho0 = 2.6867774e19;
    # Avogadro constant   
    NA = 6.02214179e23;    

class Xe:

    X0 = 1530. # xenon radiation length  * pressure in cm * bar
    mass_amu = 131.293;        # mass of xenon in amu

class NEXT:

    def __init__(self,pgas=10.,tgas=293.15):
        self.pgas = pgas # gas pressure in atm
        self.tgas = tgas # gas temperature in Kelvel
        x0 = Xe.Lr/(self.pgas/1.01325); # radiation length of xenon
        self._load()       
        return

    def _load(self):
        # Read in the stopping power and interpolate.
        xesp_tbl = np.loadtxt("xe_estopping_power_NIST.dat");
        rho = Phys.gas_rho0*(self.pgas/(self.tgas/273.15))*(Xe.mass_amu/Phys.NA);
        e_vals = xesp_tbl[:,0];
        dEdx_vals = xesp_tbl[:,1];
        e_vals = np.insert(e_vals,0,0.0);
        dEdx_vals = np.insert(dEdx_vals,0,dEdx_vals[0]);
        self.xesp = interp1d(e_vals,dEdx_vals*rho,kind='cubic')
        return 

    def stepsize(self,ene0,deltaene=0.050):
        # Determine the distance of this step.
        enef = max(0.,ene0-deltaene)
        dis = integrate.quad(lambda ee: 1./self.xesp(ee),enef,ene0,limit=1000)[0]
        return dis


def momentum(ene, mass = 0.511):
    """ return the momentum from a given energy 
    """
    p = sqrt((ene+mass)**2-mass**2)
    return p

class MS:

    @staticmethod
    def theta(p, dis, mass= 0.511):
        """ return the theta MS angle for a particle with mass, momentum (p) 
        that travel a distance (dis) in radiation units
        """
        ene = sqrt(mass**2+p**2)
        beta = p/ene
        tms = (13.6)/(p*1.*beta)*sqrt(dis*1.)*(1+0.038*log(dis*1.))
        #print ' thetams ene, beta, theta ',ene,beta,tms
        return tms

    @staticmethod
    def QMatrix(p,zdis,tx,ty,mass=0.511,X0=1.):
        """ returns the MS cov matrix in a reference frame (x,y,z) of a particle with momentum p and mass, traversing a material with X0 radiation lenght.
        Returns the Cov matrix for a vector [x,y,tx,ty] in the (x,y,z) reference frame, where tx,ty are the particle slopes respect z.        
        """
        norm = 1.+tx*tx+ty*ty
        cost = sqrt(1./norm)
        dis = zdis/cost
        theta0 = MS.theta(p,dis/X0,mass)
        theta02 = theta0*theta0
        ctx = norm*(1+tx*tx)
        cty = norm*(1+ty*ty)
        ctxy = norm*tx*ty
        Q = KFMatrixNull(4,4)
        Q[0,0] = zdis*zdis*theta02*ctx/3.
        Q[0,1] = zdis*zdis*theta02*ctxy/3.
        Q[0,2] = zdis*theta02*ctx/2.
        Q[0,3] = zdis*theta02*ctxy/2.
        Q[1,1] = zdis*zdis*theta02*cty/3.
        Q[1,2] = zdis*theta02*ctxy/2.
        Q[1,3] = zdis*theta02*cty/2.
        Q[2,2] = theta02*ctx
        Q[2,3] = theta02*ctxy
        Q[3,3] = theta02*cty
        for i in range(4):
            for j in range(i+1,4):
                Q[j,i] = Q[i,j]
        return Q
            
    @staticmethod
    def QMatrix0(p,zdis,mass=0.511,X0=1.):
        """ returns the multiple scattering Cov matriz for a particle of mass with a total p moving in the z direction a distance (zdis) and in a material of X0 radiation length
        The cov matrix is associated to [x,y,thetax,thetay],
        where x,y are the orthogonal coordinates to the z-direction
        """
        return MS.QMatrix(p,zdis,0.,0.,mass,X0)
    
    @staticmethod
    def Qrandom(Q):
        """ generates random numbers according with the cov matrix (Q)
        """
        return Random.cov(Q)

    @staticmethod
    def random(p,d,X0=1.,mass=0.511):
        """ generate multiple scattering random variables in x (transverse direction), theta 
        for a particle of mass and total momentum p that transverses a distance d
        in a medium with radiation lenght X0
        """
        theta0 = MS.theta(p,d/X0,mass)
        rho = sqrt(3.)/2.
        z1,z2 = random.gauss(0.,1.),random.gauss(0.,1.)
        y = z1*d*theta0/sqrt(12.)+z2*d*theta0/2.
        theta = z2*theta0
        #print 'ramdomms ',y,theta
        return y,theta

if __name__ == "__main__":
    
    """ check the MS random generation
    using simple method (x,theta) and Q cov method
    both should give same results
    """
    
    from alex.rootsvc import ROOTSvc 

    p = 2.5   # MeV
    X0 = 150. # mm
    dis = 5. # mm

    root = ROOTSvc(fname='ck_nextphysdat.root')
    root.h1d('x',100,-10.,10.)
    root.h1d('tx',100,-10.,10.)
    root.h2d('xtx',100,-10,10.,100,-10.,10.)

    root.h1d('x1',100,-10.,10.)
    root.h1d('tx1',100,-10.,10.)
    root.h2d('xtx1',100,-10,10.,100,-10.,10.)
    
    for i in range(1000):

        Q = MS.QMatrix0(p,dis,X0=X0)
        z = MS.Qrandom(Q)
        x,tx = MS.random(p,dis,X0=X0)
        
        root.fill('x',x)
        root.fill('tx',tx)
        root.fill('xtx',x,tx)

        x1,tx1 = z[0],z[2]
        root.fill('x1',x1)
        root.fill('tx1',tx1)
        root.fill('xtx1',x1,tx1)
    
    root.close()
        
    
