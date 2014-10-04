from alex.alex import IAlg
from alex.alex import Alex
from alex.rootsvc import ROOTSvc
import random
from KFBase import KFVector, KFMatrix, KFMatrixNull, KFMatrixUnitary
from kfilter import KFData,KFNode
from nextkfilter import KFLineFilter

"""

Alex example to check KFNextFilter

"""

class KFCheck(IAlg):
    """ Algorithm to generate SL tracks, fit them with KFNextFilter and fill histos
    """

    def define(self):
        self.radlen = 0.
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
        #xs = map(lambda i: self.tx*i*self.zdis,range(self.nhits))
        # true xs positions
        #ys = map(lambda i: self.ty*i*self.zdis,range(self.nhits)) 
        # true ys positions
        zs = map(lambda i: i*self.zdis+self.zz0,range(self.nhits)) 
        # true zs positions
        #xss = map(lambda x: x+random.gauss(0.,self.xres),xs) 
        # measured xs positions
        #yss = map(lambda x: x+random.gauss(0.,self.yres),ys) 
        # measured xs positions
        #print ' xss ',xss
        #print ' yss ',yss
        #print ' zs ',zs

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
        kf = KFLineFilter(radlen=self.radlen)
        kf.sethits(nullhits)

        x0 = KFVector([0.,0.,0.,0.,2.5])
        C0 = KFMatrixNull(5,5)
        state0 = KFData(x0,C0,0.)
    
        hits,states = kf.generate(state0)
        #print ' hits ',hits
        
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
        kf.filter(state0)  
        kf.smoother()
        #print ' nodes ',kf.nodes

        # draw results
        #-------------------

        # for a given node (ihit)
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
        for i in range(self.nhits):
            node = kf.nodes[i]
            res = node.residual('smooth')
            x,sx = node.param('smooth',0)
            y,sy = node.param('smooth',1)
            tx,stx = node.param('smooth',2)
            ty,sty = node.param('smooth',3)
            root.fill("hxi",i,x)
            root.fill("hyi",i,y)
            root.fill("htxi",i,tx)
            root.fill("htyi",i,ty)
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


if __name__ == '__main__':

    alex = Alex()
    alex.nevts = 1000
    alex.evt
    root = ROOTSvc('root','ck_anextkfilter2.root')
    alex.addsvc(root)

    alg = KFCheck('kfcheck')
    alg.radlen = 1500.
    alg.imports.append('root')

    alex.addalg(alg)

    alex.run()
    
    

