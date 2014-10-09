from alex.alex import IAlg
from alex.alex import Alex
from alex.rootsvc import ROOTSvc
import random
from KFBase import KFVector, KFMatrix, KFMatrixNull, KFMatrixUnitary
from kfilter import KFData,KFNode
from nextkfilter import KFZLineFilter,KFNextGenerator
from nextkfnoiser import NEXT,MSNoiser,ELoss
from math import *

"""

Alex example to check KFNextFilter

"""

class KFHistos(IAlg):

    def define(self):
        self.nhits = 20

    def initialize(self):
       
        self.root.h1d('knhits',self.nhits+1,0.,self.nhits+1)       
        self.root.h1d('kfchi',100,0.,5.*self.nhits)
        self.root.h1d('kschi',100,0.,5.*self.nhits)

        self.root.h2d('hix',self.nhits,0,self.nhits,100,-10.,10.)
        self.root.h2d('hiy',self.nhits,0,self.nhits,100,-10.,10.)
        self.root.h2d('hiz',self.nhits,0,self.nhits,100,-5.,20.)
        self.root.h2d('hitx',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('hity',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('hiee',self.nhits,0,self.nhits,100,0.,3.)
        self.root.h2d('hifchi',self.nhits,0,self.nhits,100,0.,10.)
        self.root.h2d('hischi',self.nhits,0,self.nhits,100,0.,10.)
        self.root.h2d('hixpul',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('hiypul',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('hitxpul',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('hitypul',self.nhits,0,self.nhits,100,-5.,5.)
        
        self.root.h1d('hx',100,-10.,10.)
        self.root.h1d('hy',100,-10.,10.)
        self.root.h1d('hz',100,-5.,20.)
        self.root.h1d('htx',100,-5.,5.)
        self.root.h1d('hty',100,-5.,5.)
        self.root.h1d('hee',100,0.,3.)
        self.root.h1d('hfchi',100,0.,10.)
        self.root.h1d('hschi',100,0.,10.)
        self.root.h1d('hxpul',100,-5.,5.)
        self.root.h1d('hypul',100,-5.,5.)
        self.root.h1d('htxpul',100,-5.,5.)
        self.root.h1d('htypul',100,-5.,5.)

    def execute(self):

        kf = self.evt['kf']
        if (not kf): return False
        
        nhits = len(kf)
        fchi = kf.chi2('filter')
        schi = kf.chi2('smooth')
        self.root.fill('knhits',nhits)
        self.root.fill('kfchi',fchi)
        self.root.fill('kschi',schi)

        for i in range(nhits):
            node = kf.nodes[i]
            state = node.getstate('true')
            rx,ry,rtx,rty,ree = state.vec
            rz = state.zrun
            node = kf.nodes[i]
            res = node.residual('smooth')
            x,sx = node.param('smooth',0)
            y,sy = node.param('smooth',1)
            tx,stx = node.param('smooth',2)
            ty,sty = node.param('smooth',3)
            ee,see = node.param('smooth',4)
            root.fill("hix",i,x)
            root.fill("hiy",i,y)
            root.fill("hiz",i,rz)
            root.fill("hitx",i,tx)
            root.fill("hity",i,ty)
            root.fill("hiee",i,ee)
            root.fill('hixpul',i,(x-rx)/sx)
            root.fill('hiypul',i,(y-ry)/sy)
            root.fill('hitxpul',i,(tx-rtx)/stx)
            root.fill('hitypul',i,(ty-rty)/sty)
            fchi = node.getchi2('filter')
            schi = node.getchi2('smooth')
            root.fill('hifchi',i,fchi)
            root.fill('hischi',i,schi)

            root.fill("hx",x)
            root.fill("hy",y)
            root.fill("hz",rz)
            root.fill("htx",tx)
            root.fill("hty",ty)
            root.fill("hee",ee)
            root.fill('hxpul',(x-rx)/sx)
            root.fill('hypul',(y-ry)/sy)
            root.fill('htxpul',(tx-rtx)/stx)
            root.fill('htypul',(ty-rty)/sty)
            fchi = node.getchi2('filter')
            schi = node.getchi2('smooth')
            root.fill('hfchi',fchi)
            root.fill('hschi',schi)
    
        return True

class KFZLine(IAlg):
    """ Algorithm to generate SL tracks, fit them with KFNextFilter and fill histos
    """

    def define(self):
        self.radlen = -1.
        self.eloss = False
        self.x0 = 0.
        self.y0 = 0.
        self.z0 = 0.
        self.tx = 0.
        self.ty = 0.
        self.E0 = 2.5
        self.nhits = 10
        self.zz0 = 1.
        self.zdis = 1.
        self.xres = 0.01
        self.yres = 0.01
        return

    def execute(self):

        # Prepare Generation
        #----------------
        # true zs positions
        zs = map(lambda i: i*self.zdis+self.zz0,range(self.nhits)) 

        # H matrix 
        V = KFMatrixNull(2) 
        V[0,0]=self.xres*self.xres
        V[1,1]=self.yres*self.yres
        #print ' V ',V

        # Generate Track
        #----------------
        #empty this
        nullhits = map(lambda z: KFData(KFVector([0,0]),V,z),zs)
        
        # create NEXT Kalman Filter
        next = NEXT()
        msnoiser = next.msnoiser
        eloss = next.eloss
        if (self.radlen>=0):
            msnoiser = MSNoiser(self.radlen)
        if (not self.eloss): eloss=None
        kf = KFZLineFilter(msnoiser,eloss)
        kf.sethits(nullhits)

        x0 = KFVector([self.x0,self.y0,self.tx,self.ty,self.E0])
        C0 = KFMatrixNull(5,5)
        state0 = KFData(x0,C0,0.)

        # generate hits and states into nodes
        knodes = kf.generate(state0)
        # set the measured nodes in the kf
        kf.setnodes(knodes)
 
        # Fit track
        #------------------

        # fit!
        C1 = KFMatrixUnitary(5)
        state0 = KFData(x0,C1,0.)
        ok,fchi,schi = kf.fit(state0)  
        if (not ok):
            print " Fit Failed!"
            return False
        print ' Fit chi2 ',fchi,schi

        # set the kalmanfilter
        self.evt['kf'] = kf
        return True

class KFNextGen(IAlg):

    def define(self):
        self.deltae = 0.01 # 10 keV
        self.nhits = 20
        self.zz0 = 0.2
        self.dz = 0.2
        self.xres = 0.01
        return

    def execute(self):

        kfgen = KFNextGenerator(self.deltae)
        state0 = [0.,0.,0.,0.,0.,1.,2.5]
        states = kfgen.generate(state0)
        
        zs = map(lambda i: self.zz0+i*self.dz,range(self.nhits))
        
        #zstates = KFNextGenerator.sample(states,zs)
        #zhits = KFNextGenerator.hits(zstates,self.xres)
        
        znodes = KFNextGenerator.kfnodes(states,zs,self.xres)

        next = NEXT()
        kffit = KFZLineFilter(next.msnoiser,next.eloss)
        kffit.setnodes(znodes)
        C0 = KFMatrixUnitary(5)
        zstate0 = KFData(KFVector([0.,0.,0.,1.,2.5]),C0,0.)
        ok,fchi,schi = kffit.fit(zstate0)

        self.evt['gen/states'] = states
        self.evt['gen/znodes'] = znodes
        if (ok):
            self.evt['kf'] = kffit

        return ok
        

if __name__ == '__main__':

    alex = Alex()
    alex.nevts = 1000
    root = ROOTSvc('root','ck_anextkfilter.root')
    alex.addsvc(root)

    nhits = 40

    #kalg = KFZLine('kfline')
    #kalg.radlen= 1500. # for radlen
    #kalg.eloss = True # for no enery loss
    #kalg.nhits = nhits
    #kalg.eloss = True
    #kalg.zz0 = 0.2
    #kalg.zdis = 0.2

    kalg = KFNextGen('kfnextgen')
    kalg.nhits = nhits
    
    halg = KFHistos('kfhistos')
    halg.imports.append('root')
    halg.nhits = nhits

    alex.addalg(kalg)
    alex.addalg(halg)

    alex.run()
    
    

