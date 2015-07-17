from ROOT import *
import numpy as np
from array import array 

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



N = 165
h = TH2D("h","h",N+10,0,N+10,100,0,.5)
dx =TH2D("dx","dx",N,0,N,1000,0,10)
#x,y,z,ux,uy,uz,ee = [np.linspace(0,1,1000)*0. for i in range(7)]
x,y,z,ux,uy,uz,ee = [array('f',[0.]*N) for i in range(7)]
ff = TFile("TestGen.root","READONLY")

tt = ff.Get("TestGen")
tt.SetBranchAddress('x',x)
tt.SetBranchAddress('y',y)
tt.SetBranchAddress('z',z)
tt.SetBranchAddress('ux',ux)
tt.SetBranchAddress('uy',uy)
tt.SetBranchAddress('uz',uz)
tt.SetBranchAddress('ee',ee)

for j in range(int(tt.GetEntries())):
    tt.GetEntry(j)
    for i in range(1,N):
#        print zip(x,y,z,ux,uy,uz,ee)[i]
#        continue
        dr = ((x[i]-x[i-1])**2+(y[i]-y[i-1])**2+(z[i]-z[i-1])**2)**0.5
        if dr:
            h.Fill(i,(ee[i-1]-ee[i])/dr)
            dx.Fill(i,dr)
    e = map(lambda g : ee[g-1]-ee[g] ,range(len(ee)))
    a = Plot4D(x,y,z,e)
    print j, -min(x)+max(x), -min(y)+ max(y), -min(z) + max(z)
    raw_input()
    del a


