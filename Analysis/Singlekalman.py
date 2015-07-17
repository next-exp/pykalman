from atrackgen import *
from akfbeta import *
from math import *
from messenger import *
import numpy as np
from Ana import *


def Plot4D( x, y, z, t, markerstyle = 20, markersize = 1 ):
    '''
        Plot a 3D dataset (x,y,z) with an extra color coordinate (t).
    '''
    data = array( 'd', [0.] * 4 )
    tree = TTree('DummyTree','DummyTree')
    tree.Branch('xyzt', data, 'x/D:y:z:t')

    for datai in zip(x,y,z,t):
        data[0], data[1], data[2], data[3] = datai
        tree.Fill()
    tree.SetMarkerStyle( markerstyle )
    tree.SetMarkerSize( markersize )
    c = TCanvas()
    tree.Draw('x:y:z:t','','zcol')
    return c, tree
    


def PasandoEstados(Traza):
    states = [x for x in Traza.states if x != (0.,0.,0.,0.,0.,1.,0.)]

    return states




class ReadStates(IAlg):
    '''
    Con esto leemos del dic evt las trazas y creamos los estados, que son listas de las variables pero puestas punto por punto
    '''
    def __init__(self,input):
    
        IAlg.__init__(self,'ReadStates')
        self.input = input


    def execute(self):
               
        states = PasandoEstados(self.evt[self.input])
        states = map( list,states)
        for state in states:
            
            for i in range(3):
                state[i] *= 0.1 

        self.evt['sim/states'] = [states]

        return True

def Voxel(track): #Traza orenada en hits
    '''
    Metodo de Brais para la voxelizacion. En algun momento cambiamos para el de gonzalo que esta en Ana.py
    '''   
    
    x,y,z,ux,uy,uz,E = zip(*track)
    
    n = len(x)
    Eaux = E[0]
    Et = []
    for i in range(n-1):
        Et.append(Eaux-E[i+1])
        Eaux = E[i+1]
    Et.append(Eaux)
    
    zmin = min(z)
    zmax = max(z)
    desp = 0
    zplanes = [i/10.+desp for i in range(10*int(round(zmin)-10)-1,10*int(round(zmax))+10)]
    z0 = z[0]
    ok = False
    aux = 0

    while not ok:
        
        if z[0]<zplanes[aux]:
            ok = True
        else:
            aux += 1
    Eantes = E[0]
    vx,vy,vz,vux,vuy,vuz,vE,vEt = [],[],[],[],[],[],[],[Eantes]
    auxx,auxy,auxE = [],[],[]
    for i in range(n):
        
        if z[i]<=zplanes[aux] and z[i]>=zplanes[aux-1]:
            auxx.append(x[i])
            auxy.append(y[i])
            auxE.append(E[i])
            
        elif z[i]>zplanes[aux]:
            aux +=  int( ( z[i] - zplanes[aux] )*10 )
            vx.append(np.mean(np.array(auxx)))
            vy.append(np.mean(np.array(auxy)))
            vz.append(zplanes[aux])
            vux.append(ux[i])
            vuy.append(uy[i])
            vuz.append(uz[i])
            vE.append(auxE[0]-E[i])
            #print 'mas',auxE,E[i]
            vEt.append(Eantes-auxE[0]+E[i])
            Eantes = Eantes-auxE[0]+E[i]
            auxx,auxy,auxE = [],[],[]
            auxx.append(x[i])
            auxy.append(y[i])
            auxE.append(E[i])
        elif z[i]<zplanes[aux-1]:
            aux += -1+int( ( z[i] - zplanes[aux] )*10 )
            vx.append(np.mean(np.array(auxx)))
            vy.append(np.mean(np.array(auxy)))
            vz.append(zplanes[aux-1])
            vux.append(ux[i])
            vuy.append(uy[i])
            vuz.append(uz[i])
            vE.append(auxE[0]-E[i]) 
            #print 'merino',auxE,E[i] 
            vEt.append(Eantes-auxE[0]+E[i]) 
            Eantes = Eantes-auxE[0]+E[i]
            auxx,auxy,auxE = [],[],[]
            auxx.append(x[i])
            auxy.append(y[i])
            auxE.append(E[i])
   
        else:
            print z[i],zplanes[aux-1]
            print 'problem'
    
    if len(auxx):
        vx.append(np.mean(np.array(auxx)))
        vy.append(np.mean(np.array(auxy)))
        vz.append(zplanes[aux])
        vux.append(ux[i])
        vuy.append(uy[i])
        vuz.append(uz[i])
        vE.append(auxE[0]-E[i]) 
        
        
    qv = Plot4D(vx,vy,vz,vEt)
    raw_input('aaaa')
    
    for i in range(1,len(vz)-1):
        if vz[i]==vz[i-1]:
            print '*******************************************************************','\n'
            print i,len(vz),vz[i],'\n'
            print z,'\n'
            print vz,'\n'
            print zplanes,'\n'
            print '*******************************************************************'
    
    return [zip(vx,vy,vz,vE)[:-3]],[zip(vx,vy,vz,vux,vuy,vuz,vEt)[:-3] ]




class cVoxel(IAlg):
    """ Algorithm to generate digits, separated a distance in z (dz)
    """

    def define(self):
        
        return

    def initialize(self):
        return True

    def execute(self):
        states = self.evt['sim/states']
        if (not states): return False
        digs,zsts = Digitalize(states[0])
        self.evt['sim/digits']=digs
        self.evt['sim/zstates']=zsts
        data1 = map(lambda dig: (len(dig),dig[0],dig[-1]),digs)
        data2 = map(lambda zst: (len(zst),zst[0],zst[-1]),zsts)
        self.msg.verbose(self.name,' digits ',data1)
        self.msg.verbose(self.name,' zstats ',data2)    
        return True


class GonzaloKalman(IAlg):
    '''
        Gonzalo's Kalman Filter cheking version 
    '''
    def define(self):
        self.opath = 'Gonza'
        return

    def initialize(self):
        return True

    def execute(self):
        states = self.evt['sim/states']  
        if (not states): return False
        x, y, z, ux, uy, uz, E = zip(*states[0])
        digits,zsts= Digitalize( states[0])
        #print zsts
        hits = Hits( [digits] )
        nodes = Nodes( hits)#,states = zsts)
        AddTotalEnergy( [nodes] )
        c, kf = Fit( [nodes] )  
        kfs = []
        kfs.append(kf)
        #map(lambda kf: (kf.cleanchi2('filter'),kf.cleanchi2('smooth')),kfs)
        self.evt[self.opath]=kfs

        return True


def ck_gen():

    alex = Alex()
    alex.nevts = 1000
    
    root = ROOTSvc('root','SimpleBetawMyVoxel_Irene_15_0T_10.root')
    alex.addsvc(root)
    alex.msg.level = Message.Crucial

    RR = lambda fname : TreeReader(fname,tname='ktree') #myone
    
    reader = Reader(name='reader',freader=RR)
    #read.fnames = ['trackgen1.root','trackgen2.root','trackgen3.root']
    reader.fnames = ['SE_15atm_0T.root']# ['DatosAMirar.root']



    leer = TrackRecover('readmyone')

    Brais = ReadStates('readmyone')
        
    alex.addalg(reader)
    alex.addalg(leer)
    alex.addalg(Brais)

    hisgenstates = HistosGenerateStates('hisgenstates')  
    hisgenstates.imports.append('root')
    alex.addalg(hisgenstates)

    # generate hits
    agendigits = cVoxel('agendigits')
    alex.addalg(agendigits)

    # create nodes
    agennodes = CreateNodes('agennodes')
    alex.addalg(agennodes)


    # fit
    #-----------------

    # fit in the forward direction
    akf = FitNodes2('akf')
    akf.path = 'rec/nodes'
    akf.opath = 'rec/kfs'
    alex.addalg(akf)

    # histogram the forward fit
    hisakf = HistosKFFit('hisakf')
    hisakf.imports.append('root')
    hisakf.prefix='kf_'
    hisakf.path='rec/kfs'
    alex.addalg(hisakf)

    # histogram the chi2-scan of the forward fit
    hischiscan = HistosKFChiScan('hischiscan')
    hischiscan.imports.append('root')
    hischiscan.path = 'rec/kfs'
    alex.addalg(hischiscan)

    # reverse the nodes to do a reverse fit
    arevnodes = ReverseNodes('arevnodes')
    arevnodes.path = 'rec/nodes'
    arevnodes.opath = 'rec/nodes/rev'
    alex.addalg(arevnodes)

    # do reverse fit
    akfrev = FitNodes('akfrev')
    akfrev.path = 'rec/nodes/rev'
    akfrev.opath = 'rec/kfs/rev'
    alex.addalg(akfrev)

    # histogram the reverse fit
    hisakfrev = HistosKFFit('hisakfrev')
    hisakfrev.imports.append('root')
    hisakfrev.prefix='kfrev_'
    hisakfrev.path='rec/kfs/rev'
    alex.addalg(hisakfrev)

    # histogram the chi2-scan of the reverse fit
    hischiscanrev = HistosKFChiScan('hischiscanrev')
    hischiscanrev.imports.append('root')
    hischiscanrev.path='rec/kfs/rev'
    hischiscanrev.prefix='chiscanrev_trk'
    alex.addalg(hischiscanrev)

    # Compare forward/reverse fit
    #-----------------------------

    # do the analysis of the forward/reverse fit
    anarev = HistosAnaReverse('anrev')
    anarev.imports.append('root')
    alex.addalg(anarev)
    
    alex.run()    

def Gon_gen():
    '''
    Not used
    '''
    alex = Alex()
    alex.nevts = 1000
    root = ROOTSvc('root','MyVoxelPerCent.root')
    alex.addsvc(root)
    alex.msg.level = Message.Fatal

    RR = lambda fname : TreeReader(fname,tname='myone')
    
    reader = Reader(name='reader',freader=RR)
    #read.fnames = ['trackgen1.root','trackgen2.root','trackgen3.root']
    reader.fnames = ['DatosAMirar.root']



    leer = TrackRecover('readmyone')

    Brais = ReadStates('readmyone')
        
    alex.addalg(reader)
    alex.addalg(leer)
    alex.addalg(Brais)

    hisgenstates = HistosGenerateStates('hisgenstates')  
    hisgenstates.imports.append('root')
    alex.addalg(hisgenstates)

    # generate hits
    Kalman = GonzaloKalman('Kalman')
    alex.addalg(Kalman)

    
    # histogram the forward fit
    hisakf = HistosKFFit('hisakf')
    hisakf.imports.append('root')
    hisakf.prefix='kf_'
    hisakf.path='Gonza'
    #alex.addalg(hisakf)

    

    alex.run()        


    
if __name__ == '__main__':

    #Gon_gen()
    ck_gen()