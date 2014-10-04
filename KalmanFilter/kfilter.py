from KFBase import KFVector, KFMatrix, KFMatrixNull, KFMatrixUnitary, Random 
from math import *

"""
Generic implemantation of a Kalman Filter

Users needs to provide a KFPropagate and KFFilter class to add KFNodes 

"""

class KFData(object):
    """ class to store an vector data and its cov matrix
    it is associated to a running parameters (zrun)
    """
    
    def __init__(self,vec,cov,zrun):
        """ create KFData with a vector a cov matrix and the zrun parameter
        """
        self.vec = KFVector(vec)
        self.cov = KFMatrix(cov)
        self.zrun = zrun
        return

    def __str__(self):
        """ convert to str a KFData object
        """
        s = ' KFData '
        s += 'vector: '+str(self.vec)+'\n'
        s += 'matrix: '+str(self.cov)+'\n'
        s += 'zrun: '+str(self.zrun)
        return s

    def __repr__(self):
        """ convert to str a KFData object
        """
        return str(self)

    def copy(self):
        """ copy a KFData object
        """
        return KFData(self.vec,self.cov,self.zrun)
    
class KFPropagate(object):    
    """ virtual class to define a system and the propagation of a state
    """

    def propagate(self,state,zrun):
        """ virtual method, propagate an state to a new running position (zrun),
        it return the propagated state, F and Q matrices
        """
        return None,None,None
        
    def stepsize(self,state):
        """ returns the zrun position for the next step starting from this state
        """
        return None

class KFNode(object):
    """ note to store a mesurements and the kalman states
    It does filter and smooth functions
    """

    names = ['none','pred','filter','smooth','rpred','rfilter']
    
    def __init__(self,hit,hmatrix):
        """ constructor with a hit (KFData) with the measurment and the HMatrix
        """
        self.hit = hit
        self.hmatrix = KFMatrix(hmatrix)
        self.zrun = hit.zrun
        self.states = {}
        self.chi2 = {}
        self.status = 'none'
        return

    def __str__(self):
        """ convert a KFNode into a string
        """
        s = 'hit '+str(self.hit)+'\n'
        s+= 'states '+str(self.states)+'\n'
        s+= 'chi2 '+str(self.chi2)
        return s

    def __repr__(self):
        """ convert a KFNode into a string
        """
        return str(self)

    def setstate(self,name,state):
        """ set an state into the node ('pred','fiter','rfilter','smooth')
        """
        if (name not in KFNode.names):
            print ' state name  ',name,' not in KNode!'
        self.states[name]=state.copy()
        self.status = name
        return

    def setchi2(self,name,chi2):
        """ set a chi2 value into the node ('pred','fiter','rfilter','smooth')
        """
        if (name not in KFNode.names):
            print ' state name  ',name,' not in KNode!'
        self.chi2[name]=chi2
        return

    def getstate(self,name):
        """ get the state with name ('pred','fiter','rfilter','smooth')
        """
        return self.states[name]

    def getchi2(self,name):
        """ get the chi2 with name ('pred','fiter','rfilter','smooth')
        """
        return self.chi2[name]

    def residual(self,name):
        """ get the residual from a state  (m - H x)
        """
        state = self.getstate(name)
        m = self.hit.vec 
        x = state.vec
        res = m - self.hmatrix*x
        return res

    def param(self,name,i):
        """ get the param value and error of index i in state with name
        ('pred','fiter','rfilter','smooth')
        """
        state = self.getstate(name)
        x,C = state.vec,state.cov
        xx,cc = x[i],sqrt(C[i,i])
        return xx,cc

    def generate(self,state):
        """ generate a hit from this state
        """
        x0 = state.vec
        C = state.cov
        sx = Random.cov(C)
        x = x0+sx
        #print ' x0 ',x0,' sx ',sx,' x ',x
        n,m = C.M.shape
        C0 = KFMatrixNull(n,m)
        gstate = KFData(x,C0,self.zrun)
        H = self.hmatrix
        m0 = H*x
        V = self.hit.cov
        sm = Random.cov(V)
        m = m0+sm
        #print ' m0 ',m0,' sm ',sm,' m ',m
        hit = KFData(m,V,self.zrun)
        #print ' initial state ',state
        #print ' generate state ',gstate
        #print ' generate hit ',hit
        return hit,gstate

    def predict(self,state):
        """ state is the propagated state to this node
        return the filter state at this node and the chi2
        """
        #print ' KFNode predict at ',self.zrun
        x = state.vec
        C = state.cov
        Ci = C.Inverse()
        #print ' C ',C
        #print ' Ci ',Ci
        m = self.hit.vec
        V = self.hit.cov
        Vi = V.Inverse()
        #print ' V ',V
        #print ' Vi ',Vi
        H = self.hmatrix   
        HT = H.Transpose()
        #print ' H ',H
        #print ' HT ',HT
        Cfi = Ci+(HT*Vi)*H
        Cf = Cfi.Inverse()
        #print ' Cfi ',Cfi
        #print ' Cf ',Cf
        xf = Cf*(Ci*x+HT*(Vi*m))
        #print ' xf',xf
        fstate = KFData(xf,Cf,self.zrun)
        mres = m-H*xf
        mrest = mres.Transpose()
        xres = xf-x
        xrest = xres.Transpose()
        #print ' mres ',mres,' xres ',xres
        chi2 = mrest*Vi*mres+xrest*Ci*xres
        chi2 = chi2[0]
        #print ' filtered state ',fstate
        #print ' chi2 ',chi2
        return fstate,chi2
        
    def smooth(self,node1):
        """ node is the next node already smoothed
        return the smother state at this node and the chi2
        """
        fstate = self.getstate('filter')
        xf = fstate.vec
        Cf = fstate.cov
        pstate1 = node1.getstate('pred')
        xp1 = pstate1.vec
        Cp1 = pstate1.cov
        Cp1i = Cp1.Inverse()
        F = self.F
        FT = F.Transpose()
        A = Cf*(FT*Cp1i)
        AT = A.Transpose()
        sstate1 = node1.getstate('smooth')
        xs1 = sstate1.vec
        Cs1 = sstate1.cov
        xs = xf + A*(xs1-xp1)
        Cs = Cf + A*(Cs1-Cp1)*AT
        sstate = KFData(xs,Cs,self.zrun)
        m = self.hit.vec
        V = self.hit.cov
        H = self.hmatrix
        HT = H.Transpose()
        res = m-H*xs
        rest = res.Transpose()
        R = V-H*(Cs*HT)
        Ri = R.Inverse()
        chi2 = rest*(Ri*res)
        chi2 = chi2[0]
        #print " smooth state ",sstate
        #print " smooth chi2 ",chi2
        return sstate,chi2
        

class KFFilter(object):
    """ Base clase for Kalman Filter
    It requires the list of nodes and a propagator (KFPropagate)
    """

    def __init__(self,nodes=[],propagator=[]):
        """ empty constructor
        nodes are the KFNodes to fit
        propagator is a KFPropagate instance to propagate states to nodes
        """
        self.nodes = nodes
        self.propagator = propagator
        self.status = 'none'
        return

    def generate(self,state0):
        state = state0.copy()
        hits = []
        states = []
        for node in self.nodes:
            zrun = node.zrun
            state,F,Q = self.propagator.propagate(state,zrun)
            hit,state = node.generate(state)
            hits.append(hit)
            states.append(state.copy())
        return hits,states

    def filter(self, state0):
        """ executes the Kalman Filter (go only) from using a seed state (state0)
        """
        state = state0.copy()
        tchi2 = 0
        for node in self.nodes:
            zrun = node.zrun
            state,F,Q = self.propagator.propagate(state,zrun)
            node.F = F
            node.Q = Q
            node.setstate('pred',state)
            fstate,fchi2 = node.predict(state)
            node.setstate('filter',fstate)
            node.setchi2('filter',fchi2)
            tchi2+=fchi2
            state = fstate.copy()
        self.status='filter'
        #print " total chi2 ",tchi2
        return tchi2

    def rfilter(self,state0):
        """ executes the Kalman Filter (go only) in reverse mode 
        using a seed state (state0)
        """
        state = state0.copy()
        tchi2 = 0
        ks = range(len(self.nodes))
        ks.revese()                   
        for k in ks:
            node = self.nodes[k]
            zrun = node.zrun
            state,F,Q = self.propagator.propagate(state,zrun)
            node.Fr = F
            node.Qr = Q
            node.setstate('rpred',state)
            fstate,fchi2 = node.predict(state)
            node.setstate('rfilter',fstate)
            node.setchi2('rfilter',fchi2)
            tchi2+=fchi2
        self.status='rfilter'
        #print " total chi2 ",tchi2
        return tchi2

    def smoother(self):
        """ executes the Kalman Smother (back) after the filter is applied
        """
        if (self.status != 'filter'):
            print 'Status is not Filter '
            return
        fstate = self.nodes[-1].getstate('filter')
        self.nodes[-1].setstate('smooth',fstate.copy())
        self.nodes[-1].setchi2('smooth',self.nodes[-1].getchi2('filter'))
        ks = range(0,len(self.nodes)-1)
        ks.reverse()
        tchi2 = 0.
        for k in ks:
            node = self.nodes[k]
            node1 = self.nodes[k+1]
            sstate,schi2 = node.smooth(node1)
            node.setstate('smooth',sstate)
            node.setchi2('smooth',schi2)
            tchi2+=schi2
        self.status='smooth'
        #print " total chi2 ",tchi2
        return tchi2

    def fit(self,state0):
        """ execute the full Kalman filter + smoother from a seed state (state0)
        """
        fchi2 = self.filter(state0)
        schi2 = self.smoother()
        return fchi2,schi2

    def clear(self):
        """ clears the nodes
        """
        self.nodes = []
        self.status = 'none'
        return
