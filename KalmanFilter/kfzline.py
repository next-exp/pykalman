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

def iszstate(state):
    """ return true if this state is a ZLine state 
    """
    ok = hasattr(state,'uz')
    debug('isstate ok',ok)
    return ok

def zstate(ustate):
    """ create a zstate from a ustate (x,y,z,ux,uy,uz,ee)
    into (x,y,tx,ty,ee) and uz in the pars of KFData
    """
    x0,y0,z0,ux,uy,uz,ee = ustate
    assert uz != 0.,"zline.zstate uz cos is null!"
    assert ee != 0.,"zline.zstate energy is null!"
    x = KFVector([x0,y0,ux/uz,uy/uz,ee])
    C = KFMatrixNull(5)
    zst = KFData(x,C,zrun=z0,pars={'uz':uz})
    debug('zsstate ',(ustate,zst))
    return zst


class KFZLine(KFModel):
    """ StraightLine propagation of a ZState
    """

    def __init__(self,noiser=None,eloss=None):
        """ constructor of a ZLine propoagator for a (x,y,tx,ty,ene) statte with
        para the forward(+1)/backward(-1) sense along z
        It uses a noiser for the MS and energy loss for ene parameters
        """
        self.noiser = noiser
        self.eloss = eloss 
        return

    def validstep(self,state,zrun):
        """ return true is this step is physically acceptable
        """
        ok = iszstate(state)
        assert ok,'kfzline.validstep not valid state'
        dz = zrun-state.zrun
        #ok = ok and (state.uz*dz>0.)
        ene = state.vec[-1]
        if (self.noiser):
            ok = ok and self.noiser.validstep(ene,dz)
        if (self.eloss):
            ok = ok and  self.eloss.validstep(ene,dz)
        if (not ok):
            warning('kfzline.validstep not valid step',(state.vec,zrun))
        debug('kfzline.validstep ok ',ok)
        return ok
        
    def QMatrix(self,x,zdis):
        """ Returns the Q matrix for a deltaz
        """ 
        x0,y0,tx,ty,p = x
        norm = 1.+tx*tx+ty*ty
        cost = sqrt(1./norm)
        dis = zdis/cost
        theta0 = self.noiser.theta(p,dis)
        theta02 = theta0*theta0
        ctx = norm*(1+tx*tx)
        cty = norm*(1+ty*ty)
        ctxy = norm*tx*ty
        Q = KFMatrixNull(5,5)
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
        Q[4,4]=1.e-32
        debug('ZLineMode.Q x,zdis,Q ',(x,zdis,Q))
        return Q

    def FMatrix(self,x,dz):
        """ return the FMatrix
        """
        n = x.Length()
        F = KFMatrixUnitary(n)
        F[0,2]=dz
        F[1,3]=dz
        if (self.eloss):
            tx,ty,ene0=x[2],x[3],x[4]
            icost = sqrt(1.+tx*tx+ty*ty)
            ds = dz*icost
            dene,ds = self.eloss.deltae(ene0,ds)
            rat = (ene0-dene)/ene0
            assert rat >= 0,'KFZLinePropoage.FMatrix no valid eloss %f'%rat
            F[4,4] = rat
        debug('KFNextPropagate.FMatrix x,dz,F ',(x,dz,F))
        return F