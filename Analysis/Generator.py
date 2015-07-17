from kfnext import NEXT,nextgenerator

from math import *
from random import Random
import shelve

random = Random()
update = 100


#N = 500
#fout = shelve.open( 'single.shelve' )
#fout = open( 'single.dat', 'w' )
#fout.write('### x y z ux uy uz E ###\n')

#for i in range(N):
#    if not i%update:
#        print '### =>',i
#    fout[str(i)] = genele()
#    track  = genele()
#    [ fout.write( ' '.join( map( str, state ) ) + '\n' ) for state in track ]
#    fout.write('\n')
#fout.close()

def GenerateSingleBeta( N = 100, filename = 'singlebeta.shelve', **options ):
    fout = shelve.open( filename )
    fout['N'] = N
    P     = options.get(    'P', 15   )
    E     = options.get(    'E', 2.5  )
    dE    = options.get(   'dE', 0.05 )
    Emin  = options.get( 'Emin', 0.05 )
    zshot = options.get( 'zshot' )
    generator = nextgenerator( NEXT(P), dE, Emin )
    
    for i in range(N):
        if not i % update: print '### => ', i
        if zshot:
            state0 = (0.,0.,0.,0.,0.,1.,E)
        else:
            theta = acos( random.uniform(-1.,1.) )
            phi   = random.uniform( 0, 2.*pi )
            ux,uy,uz = sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)
            state0 = (0.,0.,0.,ux,uy,uz,E)
        fout[str(i)] = generator.generate( state0 )
    fout.close()

def GenerateDoubleBeta( N = 100, filename = 'doublebeta.shelve', **options ):
    fout = shelve.open( filename )
    fout['N'] = N
    P     = options.get(    'P', 15   )
    E     = options.get(    'E', 2.5  )
    dE    = options.get(   'dE', 0.05 )
    Emin  = options.get( 'Emin', 0.05 )
    generator = nextgenerator( NEXT(P), dE, Emin )
    
    for i in range(N):
        if not i % update: print '### => ', i
        E1 = random.uniform( 0., E - 2 * Emin )
        E2 = E - E1
        tracks = []
        for e in [E1,E2]:
            theta = acos( random.uniform(-1.,1.) )
            phi   = random.uniform( 0, 2.*pi )
            ux,uy,uz = sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)
            state0 = (0.,0.,0.,ux,uy,uz,e)
            tracks.append( generator.generate( state0 ) )
        fout[str(i)] = tracks
    fout.close()


GenerateSingleBeta( 100, 'sbeta.shelve', dE = 0.001, Emin = 0.000 )
#GenerateDoubleBeta( 1000 )

