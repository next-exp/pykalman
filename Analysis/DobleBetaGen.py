#from atrackgen import *
from akfbeta import *
from numpy import array as narray
from array import array


def Plot4D( x, y, z, t, markerstyle = 20, markersize = 1 ):
    data = narray( [0.] * 4 )
    tree = TTree('DummyTree','DummyTree')
    tree.Branch('xyzt', data, 'x/D:y:z:t')
    
    for datai in zip(x,y,z,t):
        data[0], data[1], data[2], data[3] = datai
        #print data
        tree.Fill()
    tree.SetMarkerStyle( markerstyle )
    tree.SetMarkerSize( markersize )
    c = TCanvas()
    tree.Draw('x:y:z:t','','zcol')
    return c, tree



class SaveDoubleBeta(IAlg):

    def define(self):
        self.E0 = 2.5 # MeV
        self.path = self.name
        return
        
    def initialize(self):
        """ book the ntuple and put it in ROOTSvc
        """
        N = 492
        labels = [('x0','F',N),('y0','F',N),('z0','F',N),
                  ('ux0','F',N),('uy0','F',N),('uz0','F',N),
                  ('ee0','F',N),('x1','F',N),('y1','F',N),('z1','F',N),
                  ('ux1','F',N),('uy1','F',N),('uz1','F',N),
                  ('ee1','F',N)]
        tree = ttree(self.name,labels)
        self.root.put(tree)
        return True
    
    def pad(self,x):
        l = len(x)
        [ x.append(0.) for i in range(100-l) ]
        return
        
    def execute(self):
        states = genbeta(self.E0)
        self.evt['sim/states'] = states
        data = map(lambda seg: (len(seg),seg[0],seg[-1]),states)
        #self.msg.info(self.name,' states ',data)
        self.evt[self.name]=[states,]
        
        x0,y0,z0,ux0,uy0,uz0,ee0 = zip(*states[0])
        x1,y1,z1,ux1,uy1,uz1,ee1 = zip(*states[1])

        tup = {'x0':x0,'y0':y0,'z0':z0,'ux0':ux0,'uy0':uy0,'uz0':uz0,'ee0':ee0,
               'x1':x1,'y1':y1,'z1':z1,'ux1':ux1,'uy1':uy1,'uz1':uz1,'ee1':ee1}
        self.root.fill(self.name,tup)
        return True
        
        
#########################




alex = Alex('alex')
alex.nevts = 200

root = ROOTSvc('root',fname='DobleBetaFino_Rand.root')
alex.addsvc(root) 



gen = SaveDoubleBeta('DobBet2')
gen.imports.append('root')#from atrackgen import *
from akfbeta import *


class SaveDoubleBeta(IAlg):

    def define(self):
        self.E0 = 2.5 # MeV
        self.path = self.name
        return
        
    def initialize(self):
        """ book the ntuple and put it in ROOTSvc
        """
        N = 492
        labels = [('x0','F',N),('y0','F',N),('z0','F',N),
                  ('ux0','F',N),('uy0','F',N),('uz0','F',N),
                  ('ee0','F',N),('x1','F',N),('y1','F',N),('z1','F',N),
                  ('ux1','F',N),('uy1','F',N),('uz1','F',N),
                  ('ee1','F',N)]
        tree = ttree(self.name,labels)
        self.root.put(tree)
        return True
    
    def pad(self,x):
        l = len(x)
        [ x.append(0.) for i in range(100-l) ]
        return
        
    def execute(self):
        states = genbeta(self.E0)
        self.evt['sim/states'] = states
        data = map(lambda seg: (len(seg),seg[0],seg[-1]),states)
        #self.msg.info(self.name,' states ',data)
        self.evt[self.name]=[states,]
        states[0].extend( [(0.,0.,0.,0.,0.,0.,0.)]*(492-len(states[0])) )
        states[1].extend( [(0.,0.,0.,0.,0.,0.,0.)]*(492-len(states[1])) )
        
        x0,y0,z0,ux0,uy0,uz0,ee0 = zip(*states[0])
        x1,y1,z1,ux1,uy1,uz1,ee1 = zip(*states[1])

        tup = {'x0':x0,'y0':y0,'z0':z0,'ux0':ux0,'uy0':uy0,'uz0':uz0,'ee0':ee0,
               'x1':x1,'y1':y1,'z1':z1,'ux1':ux1,'uy1':uy1,'uz1':uz1,'ee1':ee1}
        self.root.fill(self.name,tup)
        return True
        
        
#########################




alex = Alex('alex')
alex.nevts = 1000

root = ROOTSvc('root',fname='DobleBetaFino_Rand.root')
alex.addsvc(root) 



gen = SaveDoubleBeta('DobBet2')
gen.imports.append('root')


alex.addalg(gen)
alex.run()








alex.addalg(gen)
alex.run()





