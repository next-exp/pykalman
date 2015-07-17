
"""

Physics data and formulae for NEXT

"""
from KFBase import KFMatrix, KFMatrixNull,KFMatrixUnitary
from kffilter import KFFilter
from kfgenerator import KFGenerator
from kfzline import KFZLine
from kfnoise import MSNoise, ELoss, HPXeELoss
from math import *
import random
import numpy as np

DEBUG = False
WARNING = True

def debug(comment,arg=''):
    if (DEBUG): print comment,arg

def warning(comment,arg=''):
    if (WARNING): print comment,arg


class NEXT:
    """ NEXT helper class, computes X0 and rho for a given pressure
    create MSNoise and HPXeEloss classes
    """
    rho0_gas = 2.6867774e19; # density of ideal gas at T=0C, P=1 atm in cm^(-3)
    NA = 6.02214179e23;    # Avogadro constant   
    XeX0 = 1530. # xenon radiation length  * pressure in cm * bar
    XeX0 = XeX0
    Xemass_amu = 131.293;        # mass of xenon in amu

    #pgas = 10. # atmos
    #tgas = 295.15 # Kelvin    
    def __init__(self,pgas=15.,tgas=295.15):
        """ constructor of next for a preasure and temperature
        it creates a MS with X0 and an ELoss for electrons
        """
        self.pgas = pgas
        self.tgas = tgas
        self.X0 = NEXT.XeX0/(pgas/1.01325); # radiation length of xenon
        self.rho = NEXT.rho0_gas*(pgas/(tgas/273.15))*(NEXT.Xemass_amu/NEXT.NA);
        self.msnoiser = MSNoise(self.X0)
        self.eloss = HPXeELoss(self.rho)
        debug('NEXT pgas, X0, rho ',(pgas,self.X0,self.rho))

#-------------------------------------------------------

# basic projection matrix
H0 = KFMatrixNull(2,5); H0[0,0]=1; H0[1,1]=1
#H0 = KFMatrixNull(3,5); H0[0,0]=1.; H0[1,1]=1.; H0[2,4]=1.;

# basic measurement variace matrix
def V0(xres,eres=0.001): 
    """ creates a basic resolution matrix for mesurements
    """
    #uu = KFMatrixUnitary(3)
    uu = KFMatrixUnitary(2) 
    uu[0,0]=xres*xres; uu[1,1]=xres*xres; 
    #uu[2,2]=eres*eres
    return uu

# next parameters and ms and eloss helper methods
next0 = NEXT()

def nextgenerator(next=next0,deltae=0.005,emin=0.05):
    """ creates a generators of particle states with a given next configuration
    default NEXT 10 atms
    """
    kf = KFGenerator(next.msnoiser,next.eloss,deltae=deltae,emin=emin)
    return kf

def nextfilter(next=next0):
    """ creates a Kalman Filter for ZLine states with next configuration
    default NEXT 10 atms
    """
    zmodel = KFZLine(next.msnoiser,next.eloss)
    kf = KFFilter(model=zmodel)
    return kf

def simplegenerator(radlen,deltae=0.02,emin=0.05):
    """ creates a generator of particle states with eloss HPXe 10 atms 
    and radlen
    """
    noiser = MSNoise(radlen)
    eloss = next0.eloss
    #eloss = None
    kf = KFGenerator(noiser,eloss,deltae=deltae,emin=emin)
    return kf
    
def simplefilter(radlen):
    """ creates a Kalman Filter for ZLine states with eloss HPXe 10 atms 
    and radlen
    """
    noiser = MSNoise(radlen)
    eloss = next0.eloss
    #eloss = None
    zmodel = KFZLine(noiser,eloss)
    kf = KFFilter(model=zmodel)
    return kf

