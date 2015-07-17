from akfbeta import createzdigits as Digits, createhits as Hits, createnodes as Nodes, preparenodes as AddTotalEnergy, fitnodes as Fit
from kfnext import NEXT, nextfilter as KF
from Plotting import *
from ROOT import *
from math import *
from copy import copy, deepcopy
import shelve

def mean( x ):
    return sum(x)/len(x)

def Unify( track1, track2 ):
    return track1[::-1] + track2

def Digitalize( states, width = 0.1, xyres = 0.1 ):
    x, y, z, ux, uy, uz, E = zip(*states)
    zmin, zmax = min(z), max(z)
    if zmin == zmax:
        print 'WARNING! zmin == zmax!!!'
    xd, yd, zd, uxd, uyd, uzd, Ed = [], [], [], [], [], [], []
    z0 = zmin
    while z0 < zmax:
        z1 = z0 + width
        states_in_bin = filter( lambda state: z0 <= state[2] < z1, states )
        x_in, y_in, z_in, ux_in, uy_in, uz_in, E_in = zip(*states_in_bin)
        xd.append( mean(x_in) )
        yd.append( mean(y_in) )
        zd.append( z0 + 0.5 * width )
        uxd.append( states_in_bin[-1][3] )
        uyd.append( states_in_bin[-1][4] )
        uzd.append( states_in_bin[-1][5] )
        Ed.append( max(E_in) - min(E_in) )
        z0 = z1
    return zip(xd,yd,zd,uxd,uyd,uzd,Ed)

def Digitalize( states, width = 0.1, xyres = 0.1 ):
    states = deepcopy(states)
    
    zs = [ state[2] for state in states ]
    z0 = z00 = min(zs)

    z1 = max(zs)
    zbins = []
    while z0 <= z1:
        zbins.append( z0 )
        z0 += width
    
    nbins = len(zbins)
    def FindBin( state ):
        z = state[2]
        z0 = z00
        for i in range(nbins-1):
            if z < zbins[i+1]:
                return i
        return nbins

    xd, yd, zd, uxd, uyd, uzd, Ed = [], [], [], [], [], [], []
    Eaux = states[0][6]
    
    #print 'caca',states,'\n'
    #print '**************************'
    #print Eaux
    Etd = []
    while states:
        bin0 = FindBin( states[0] )
        upper = 1
        for i in range(1,len(states)):
            if not FindBin(states[i]) == bin0:
                break
            upper += 1
        states_in_bin = states[:upper]
        # states_in_bin,'*********************'
        states        = states[upper:]
        x_in, y_in, z_in, ux_in, uy_in, uz_in, E_in = zip(*states_in_bin)
        
        Edep = max(E_in) - states[0][-1] if len(states) else max(E_in)
        xd.append( mean(x_in) )
        yd.append( mean(y_in) )
        zd.append( mean(z_in))#zbins[upper] - 0.5 * width )
        uxd.append( states_in_bin[-1][3] )
        uyd.append( states_in_bin[-1][4] )
        uzd.append( states_in_bin[-1][5] )
        Ed.append( Edep )
#        Etd.append(Eaux-max(E_in) + min(E_in))
        Etd.append(Eaux)
        
        Eaux -= Edep

    return [zip( xd, yd, zd, Ed )],[zip( xd, yd, zd, uxd, uyd, uzd, Etd )]

def Invert(lis):
    m = max(lis)
    return [m-i for i in lis]
'''
fin = shelve.open( 'sbeta.shelve' )
N = fin['N']
for i in range(N):
    states = fin[str(i)]
    nstates = len(states)
    x, y, z, ux, uy, uz, E = zip(*states)
    #xd, yd, zd, uxd, uyd, uzd, Ed = zip(*Digitalize(states))
    dx = [ x[i] - x[i-1] for i in range(1,nstates) ]
    dy = [ y[i] - y[i-1] for i in range(1,nstates) ]
    dz = [ z[i] - z[i-1] for i in range(1,nstates) ]
    dE = [ E[i] - E[i-1] for i in range(1,nstates) ]
    dr = [ dxi*dxi + dyi*dyi + dzi * dzi for dxi,dyi,dzi in zip(dx,dy,dz) ]
    dr = map( sqrt, dr )
    
    dEdx = [ -dEi/dri for dEi, dri in zip(dE,dr) ]
    plot = Plot4D( x, y, z, Invert(E) )
    #plot2 = Plot4D( xd, yd, zd, Ed )
    C = TCanvas()
    h = Graph( range(1,nstates), dEdx, 'step', 'dE/dx (MeV/cm)', 'Bragg', 20 )
    h.SetMinimum(0)
    h.Draw('AP')
    raw_input()

'''

'''
fin = shelve.open( 'doublebeta.shelve' )
N = fin['N']
for i in range(N):
    track1, track2 = fin[str(i)]
    nstates1, nstates2 = len(track1), len(track2)
    bbtrack = Unify( track1, track2 )
    x, y, z, ux, uy, uz, E = zip(*bbtrack)
    print track1[0][-1],track2[0][-1]
    plot = Plot4D( x, y, z, Invert(E) )
    raw_input()
    del plot
'''

'''
fin = shelve.open( 'sbeta.shelve' )
N   = fin['N']

for i in range(N):
    states = fin[str(i)]
    x, y, z, ux, uy, uz, E = zip(*states)
    digits, zs = Digits( states, z )
    hits = Hits( digits )
    print digits
    nodes = Nodes( hits )
    AddTotalEnergy
    
    ( nodes )
    c, k = Fit( nodes )
    #a = PlotTrack( x, y, z, Invert(E), k, 'filter' )
    print 'Fitted'
    raw_input('Enter')
'''


