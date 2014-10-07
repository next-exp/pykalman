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

class NEXT:
    """ NEXT parameters (depends on pressure and temperature)
    """
    # density of ideal gas at T=0C, P=1 atm in cm^(-3)
    rho0_gas = 2.6867774e19;
    # Avogadro constant   
    NA = 6.02214179e23;    
    XeX0 = 1530. # xenon radiation length  * pressure in cm * bar
    Xemass_amu = 131.293;        # mass of xenon in amu
    
    def __init__(self,pgas=10.,tgas=293.15):
        self.pgas = pgas # gas pressure in atm
        self.tgas = tgas # gas temperature in Kelvel
        self.X0 = NEXT.XeX0/(self.pgas/1.01325); # radiation length of xenon
        self.rho = NEXT.rho0_gas*(self.pgas/(self.tgas/273.15))*(NEXT.Xemass_amu/NEXT.NA);

        self.eloss = ELoss(self.rho)
        self.ms = MS(self.X0)
        
        return

class ELoss:
    """ Energy loss for electrons in NEXT (depend on the density rho)
    """

    def __init__(self,rho):
        """ construction of Energy loss object
        rho is the Xe density
        """
        self.rho = rho
        # Read in the stopping power and interpolate.
        xesp_tbl = np.loadtxt("xe_estopping_power_NIST.dat");
        self.evals = xesp_tbl[:,0]
        self.dEdxvals = xesp_tbl[:,1];
        self.evals = np.insert(self.evals,0,0.0);
        self.dEdxvals = np.insert(self.dEdxvals,0,self.dEdxvals[0]);
        self.xesp = sc_interpol(self.evals,self.dEdxvals*self.rho,kind='cubic')
        # this is the dE/dx(E) function
        return 

    def isvalid(self,ene0,ds):
        ok = True
        ds = abs(ds)
        de,dz = self.deltaE(ene0,ds)
        if (dz+0.001<ds): ok = False
        if (ene0-de<=0.): ok= False
        if (DEBUG or not ok):
            print "ELoss.isvalid dz,ds,ene0,de,ok",dz,ds,ene0,de,ok
        return ok
        
    def deltaE(self,ene0,deltax=0.5):
        # temporal
        de = 0.03556*(deltax/0.5)
        de = min(ene0,de)
        deltax = (de/0.03556)*0.5
        #ff = lambda dene : deltax - self.deltax(ene0,dene)[1]
        #if (ff(0.)*ff(ene0)>0.): 
        #    return ene0,self.deltax(ene0,ene0)
        #de = sc_root(ff,0.,ene0) # find the root
        if (DEBUG): print 'Eloss.deltaE de,dx ',de,deltax 
        return de,deltax

    def deltax(self,ene0,deltaene=0.050):
        # Determine the distance of this step.
        enef = max(0.,ene0-deltaene)
        de = ene0-enef
        dis = 0.5*(de/0.03556)
        #dis = sc_integrate(lambda ee: 1./self.xesp(ee),enef,ene0,limit=1000)[0]
        if (DEBUG): print 'Eloss.deltax de,dx ',de,dis 
        return de,dis

class MS:
    """ Multiple scattering, depends on the radiation length and particle mass for a given energy
    """

    thetamax = 10.

    def __init__(self,X0,mass=0.511):
        self.X0 = X0
        self.mass = mass
        return

    def isvalid(self,p,dis):
        tms = self.theta(p,dis)
        ok = (tms<=MS.thetamax)
        if (DEBUG or not ok):
            print "MS.isvalid p,dis,tms",p,dis,tms
        return ok

    def theta(self,p,dis):
        """ return the theta MS angle for a particle with mass, momentum (p) 
        that travel a distance (dis) in radiation units
        """
        dis = abs(dis)
        if (dis<=0.): return 0.
        udis =dis/self.X0
        ene = sqrt(self.mass**2+p**2)
        beta = p/ene
        if (DEBUG): print 'MS.theta p, beta, udis, theta ',p,beta,udis
        tms = (13.6)/(p*1.*beta)*sqrt(udis*1.)*(1+0.038*log(udis*1.))
        if (DEBUG): print 'MS.theta p, beta, udis, theta ',p,beta,udis,tms
        return tms

    def Q0Matrix(self,p,dis):
        """ returns the multiple scattering Cov matriz for a particle of mass with a total p moving in the z direction a distance (zdis) and in a material of X0 radiation length
        The cov matrix is associated to [x,y,thetax,thetay],
        where x,y are the orthogonal coordinates to the z-direction
        """
        return self.QMatrix(p,dis,0.,0.)

    def QMatrix(self,p,zdis,tx,ty):
        """ returns the MS cov matrix in a reference frame (x,y,z) of a particle with momentum p and mass, traversing a material with X0 radiation lenght.
        Returns the Cov matrix for a vector [x,y,tx,ty] in the (x,y,z) reference frame, where tx,ty are the particle slopes respect z.        
        """
        norm = 1.+tx*tx+ty*ty
        cost = sqrt(1./norm)
        dis = zdis/cost
        theta0 = self.theta(p,dis)
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
        if (DEBUG): print 'MS.Q p,zdis,tx,ty,Q ',p,zdis,tx,ty,Q
        return Q
            
    
    def random(self,p,dis):
        """ generate multiple scattering random variables in x (transverse direction), theta 
        for a particle of mass and total momentum p that transverses a distance d
        in a medium with radiation lenght X0
        """
        theta0 = MS.theta(p,dis)
        rho = sqrt(3.)/2.
        z1,z2 = random.gauss(0.,1.),random.gauss(0.,1.)
        y = z1*dis*theta0/sqrt(12.)+z2*dis*theta0/2.
        theta = z2*theta0
        if (DEBUG): print 'MS.random p,dis,x,theta ',p,dis,y,theta
        return y,theta

    @staticmethod
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
        #tx = ux/uz; ty = uy/uz
        #if (not phi): phi = random.uniform(0.,2.*pi)
        #txp = tan(phi)
        #typ = (-1-txp*tx)/ty
        #vv = sqrt(1.+txp*txp+typ*typ)
        #vx = txp/vv; vy = typ/vv; vz = 1./vv
        #wx = (uy*vz-uz*vy); wy = -1.*(ux*vz-uz*vx); wz = (ux*vy-uy*vx)
        # UMatrix convert track ref system to global system
        U = KFMatrix( [[vx,wx,ux ],[vy,wy,uy],[vz,wz,uz]])
        if (DEBUG):  print 'MS.UMatrix udir, U ',udir,U
        return U        

    def XUrandom(self,p,dis,x0,udir):
        """ return a position and direction (in the global system) 
        after a random MS of a particle with momentum p
        that traverses a distance dis with a direction udir
        """
        udir.Unit()
        U = MS.UMatrix(udir)
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
            print 'MS.XYrandom p,dis,x0,udir ',p,dis,x0,udir
            print 'MS.XYrandom xt,ut ',xt,ut
            print 'MS.XYrandom dx,dxt,uf ',x0+dis*udir,U*xt,xf,uf
        return xf,uf

#--------------------------

def ck_ELoss():

    next = NEXT()
    eloss = ELoss(next.rho)
    E0 = 2.5
    dene = 0.05
    EE = E0
    xx = 0.
    print ' Delta-Ene'
    while (EE>0.0):
        dx = eloss.deltax(EE,dene)
        EE = max(EE-dene,0.)
        xx+=dx
        print ' E, x',EE,xx

    E0 = 2.5
    EE = E0
    dx = 5.
    xx = 0.
    print ' Delta-X'
    while (EE>0.0):
        dene,dx = eloss.deltaE(EE,dx)
        EE -= dene
        xx += dx
        print ' x, E ',xx,EE

    return
    
def ck_MS():
    
    """ check the MS random generation
    using simple method (x,theta) and Q cov method
    both should give same results
    """
    
    from alex.rootsvc import ROOTSvc 

    p = 2.4   # MeV
    X0 = 150. # mm
    dis = 5. # mm

    ms = MS(X0)

    root = ROOTSvc(fname='ck_nextphysdat.root')
    root.h1d('x',100,-10.,10.)
    root.h1d('tx',100,-10.,10.)
    root.h2d('xtx',100,-10,10.,100,-10.,10.)

    root.h1d('x1',100,-10.,10.)
    root.h1d('tx1',100,-10.,10.)
    root.h2d('xtx1',100,-10,10.,100,-10.,10.)
    
    for i in range(1000):

        Q = ms.QMatrix0(p,dis)
        z = Random.cov(Q)
        x,tx = ms.random(p,dis)
        
        root.fill('x',x)
        root.fill('tx',tx)
        root.fill('xtx',x,tx)

        x1,tx1 = z[0],z[2]
        root.fill('x1',x1)
        root.fill('tx1',tx1)
        root.fill('xtx1',x1,tx1)
    
    root.close()
        
def ck_XU():

    p=2.5
    dis = 10.
    ms = MS(X0=1500.)
    udir = KFVector([0.,1.,1.])
    udir.Unit()
    x0 = KFVector([0.,0.,0.])
    U = ms.XUrandom(p,dis,x0,udir)

    print 'U ',U

#-------------------    

if __name__ == "__main__":

    #ck_ELoss()
    #ck_MS()
    ck_XU()
