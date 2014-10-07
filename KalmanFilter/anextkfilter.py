from alex.alex import IAlg
from alex.alex import Alex
from alex.rootsvc import ROOTSvc
import random
from KFBase import KFVector, KFMatrix, KFMatrixNull, KFMatrixUnitary
from kfilter import KFData,KFNode
from nextkfilter import KFLineFilter,KFNextGenerator
from math import *

"""

Alex example to check KFNextFilter

"""

class KFCheck(IAlg):
    """ Algorithm to generate SL tracks, fit them with KFNextFilter and fill histos
    """

    def define(self):
        self.radlen = 0.
        self.eloss = False
        self.x0 = 0.
        self.y0 = 0.
        self.z0 = 0.
        self.tx = 0.
        self.ty = 0.
        self.nhits = 10
        self.zz0 = 1.
        self.zdis = 1.
        self.xres = 0.01
        self.yres = 0.01
        self.ihit = 4
        return

    def initialize(self):

        self.root.h1d('hnhits',50,0.,50.)
        self.root.h1d('hixpul',100,-5.,5.)
        self.root.h1d('hiypul',100,-5.,5.)
        self.root.h1d('hitxpul',100,-5.,5.)
        self.root.h1d('hitypul',100,-5.,5.)
        self.root.h1d('hifchi',100,0.,10.)
        self.root.h1d('hischi',100,0.,10.)


        self.root.h2d('hxi',self.nhits,0,self.nhits,100,-10.,10.)
        self.root.h2d('hyi',self.nhits,0,self.nhits,100,-10.,10.)
        self.root.h2d('htxi',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('htyi',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('heei',self.nhits,0,self.nhits,100,0.,3.)

        self.root.h2d('hxpul',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('hypul',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('htxpul',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('htypul',self.nhits,0,self.nhits,100,-5.,5.)
        self.root.h2d('hfchi',self.nhits,0,self.nhits,100,0.,10.)
        self.root.h2d('hschi',self.nhits,0,self.nhits,100,0.,10.)
        return True

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
        kf = KFLineFilter(radlen=self.radlen,eloss=self.eloss)
        kf.sethits(nullhits)

        x0 = KFVector([self.x0,self.y0,self.tx,self.ty,2.5])
        C0 = KFMatrixNull(5,5)
        state0 = KFData(x0,C0,0.)
    
        hits,states = kf.generate(state0)
        
        # get the true states parameters
        xs = map(lambda st: st.vec[0],states)
        ys = map(lambda st: st.vec[1],states)
        txs = map(lambda st: st.vec[2],states)
        tys = map(lambda st: st.vec[3],states)
        ees = map(lambda st: st.vec[4],states)
 
        # Fit track
        #------------------
        # set the measured hits in the kf
        kf.sethits(hits)

        # fit!
        C1 = KFMatrixUnitary(5)
        state0 = KFData(x0,C1,0.)
        ok,fchi2,schi = kf.fit(state0)  
        if (not ok):
            print " Fit Failed!"
            return False
        #print ' nodes ',kf.nodes

        # draw results
        #-------------------
        nhits = len(hits)
        self.root.fill("hnhits",nhits)

        # for a given node (ihit)
        if (self.ihit<nhits):
            node = kf.nodes[self.ihit]
            x,sx = node.param('smooth',0)
            y,sy = node.param('smooth',1)
            tx,stx = node.param('smooth',2)
            ty,sty = node.param('smooth',3)
            root.fill('hixpul',(x-xs[self.ihit])/sx)
            root.fill('hiypul',(y-ys[self.ihit])/sy)
            root.fill('hitxpul',(tx-txs[self.ihit])/stx)
            root.fill('hitypul',(ty-tys[self.ihit])/sty)
            fchi = node.getchi2('filter')
            schi = node.getchi2('smooth')
            root.fill('hifchi',fchi)
            root.fill('hischi',schi)

        # vs all nodes
        for i in range(nhits):
            node = kf.nodes[i]
            res = node.residual('smooth')
            x,sx = node.param('smooth',0)
            y,sy = node.param('smooth',1)
            tx,stx = node.param('smooth',2)
            ty,sty = node.param('smooth',3)
            ee,see = node.param('smooth',4)
            root.fill("hxi",i,x)
            root.fill("hyi",i,y)
            root.fill("htxi",i,tx)
            root.fill("htyi",i,ty)
            root.fill("heei",i,ee)
            root.fill('hxpul',i,(x-xs[i])/sx)
            root.fill('hypul',i,(y-ys[i])/sy)
            root.fill('htxpul',i,(tx-txs[i])/stx)
            root.fill('htypul',i,(ty-tys[i])/sty)
            fchi = node.getchi2('filter')
            schi = node.getchi2('smooth')
            root.fill('hfchi',i,fchi)
            root.fill('hschi',i,schi)
            #print i,sx,sy,stx,sty
        return True


class KFGen(IAlg):

    def initialize(self):

        self.root.h2d('xi',50,0.,50.,100,-20.,20.)
        self.root.h2d('yi',50,0.,50.,100,-20.,20.)
        self.root.h2d('txi',50,0.,50.,100,-10.,10.)
        self.root.h2d('tyi',50,0.,50.,100,-10.,10.)
        self.root.h2d('ei',50,0.,50.,100,0.,3.)
        self.root.h2d('zi',50,0.,50.,100,-20.,20.)

        self.root.h2d('sxi',50,0.,50.,100,0.,5.)
        self.root.h2d('stxi',50,0.,50.,100,0.,8.)


        self.root.h2d('xz',50,-10.,30.,100,-20.,20.)
        self.root.h2d('yz',50,-10.,30.,100,-20.,20.)
        self.root.h2d('txz',50,-10.,30.,100,-20.,20.)
        self.root.h2d('tyz',50,-10.,30.,100,-20.,20.)
        self.root.h2d('ez',50,-10.,30.,100,0.,3.)

        self.root.h2d('sxe',50,0.,3.,100,0.,5.)
        self.root.h2d('stxe',50,0.,3.,100,0.,8.)

    def execute(self):

        x0 = KFVector([0.,0.,0.,0.,2.5])
        E0 = 2.5; z0 = 0.
        x0 = KFVector([0.,0.,0.,0.,E0])
        C0 = KFMatrixNull(5,5)
        V = KFMatrixNull(3,3)
        xres = 0.1; eres=0.01
        V[0,0]=xres*xres; V[1,1]=xres*xres; V[2,2]=eres*eres
        H = KFMatrixNull(3,5)
        H[0,0]=1.;H[1,1]=1.;H[2,4]=1.;
        state0 = KFData(x0,C0,z0)
        
        kfgen = KFNextGenerator()
        
        hits,states = kfgen.generate(state0,H,V)

        for i in range(len(states)):
            state = states[i]
            x,y,tx,ty,ee = state.vec            
            Q = state.cov
            z = state.zrun
            self.root.fill('xi',i,x)
            self.root.fill('yi',i,y)
            self.root.fill('zi',i,z)
            self.root.fill('txi',i,tx)
            self.root.fill('tyi',i,ty)
            self.root.fill('ei',i,ee)

            self.root.fill('xz',z,x)
            self.root.fill('yz',z,y)
            self.root.fill('txz',z,tx)
            self.root.fill('tyz',z,ty)
            self.root.fill('ez',z,ee)
        
            self.root.fill('sxi',i,sqrt(Q[0,0]))
            self.root.fill('stxi',i,sqrt(Q[2,2]))

            self.root.fill('sxe',ee,sqrt(Q[0,0]))
            self.root.fill('stxe',ee,sqrt(Q[2,2]))
        return True

class KFNext(IAlg):

    def execute(self):
        
        # generation
        E0 = 2.5; z0 = 0.
        x0 = KFVector([0.,0.,0.,0.,E0])
        C0 = KFMatrixNull(5,5)
        V = KFMatrixNull(2,2)
        xres = 0.1; eres=0.01
        V[0,0]=xres*xres; V[1,1]=xres*xres; #V[2,2]=eres*eres
        H = KFMatrixNull(2,5)
        H[0,0]=1.;H[1,1]=1.;#H[2,4]=1.;
        state0 = KFData(x0,C0,z0)
        kfgen = KFNextGenerator()
        hits,states = kfgen.generate(state0,H,V)

        kffit = KFLineFilter(eloss=True)
        kffit.sethits(hits[:20])
        
        x0 = KFVector([0.,0.,0.,0.,E0])
        C0 = KFMatrixUnitary(5)
        state0 = KFData(x0,C0,z0-0.1)
        ok,fchi,schi = kffit.fit(state0)
        if (ok):
            print "Nodes 0 ",kffit.nodes[0]
            print "State 0",states[0]
            print "Nodes 1 ",kffit.nodes[5]
            print "State 5 ",states[5]

        return ok

        

if __name__ == '__main__':

    alex = Alex()
    alex.nevts = 4
    root = ROOTSvc('root','ck_anextkfilter.root')
    alex.addsvc(root)

    #alg = KFCheck('kfcheck')
    #alg.radlen = 1000.
    #alg.eloss = True
    #alg.nhits = 40
    #alg.zz0 = 0.2
    #alg.zdis = 0.2
    #alg.tx = 0.
    #alg.ty = 0.
    #alg.imports.append('root')

    #alg = KFGen('kfgen')
    #alg.imports.append('root')

    alg = KFNext('kfnext')
    alg.imports.append('root')

    alex.addalg(alg)

    alex.run()
    
    

