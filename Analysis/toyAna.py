print 'Defining environment...'

from sys import argv
from Utilities import *
from Plotting import *
from ROOT import *
from math import *

basedir  = '/Users/Gonzalo/Dropbox/Gonzalo/Facultade/NEXT/KalmanFilter'
filedir  = basedir + '/output/tracks/singlebeta'
filename = filedir + '/singlebeta_*.dat'

try:
    Nfiles = int(argv[1])
except IndexError:
    Nfiles = GetNFiles( filedir )

pi5 = 0.5 * pi
pi2 = 2.0 * pi

print 'booking...'
hx    = H1 ( 100, - 20,  20, 'hx'   , 'x (cm)'         , 'Entries'         , 'hx'     )
hy    = H1 ( 100, - 20,  20, 'hy'   , 'y (cm)'         , 'Entries'         , 'hy'     )
hz    = H1 ( 100, - 20,  20, 'hz'   , 'z (cm)'         , 'Entries'         , 'hz'     )
hux   = H1 ( 100, -  1,   1, 'hux'  , 'ux'             , 'Entries'         , 'hux'    )
huy   = H1 ( 100, -  1,   1, 'huy'  , 'uy'             , 'Entries'         , 'huy'    )
huz   = H1 ( 100, -  1,   1, 'huz'  , 'uz'             , 'Entries'         , 'huz'    )
hE    = H1 ( 100,    0, 2.5, 'hE'   , 'E (MeV)'        , 'Entries'         , 'hE'     )
hdE   = H1 ( 100,    0, 0.1, 'hdE'  , 'dE (MeV)'       , 'Entries'         , 'hdE'    )
hds   = H1 ( 100,    0,   1, 'hds'  , 'ds (cm)'        , 'Entries'         , 'hds'    )
hdx   = H1 ( 100, -  1,   1, 'hdx'  , 'dx (cm)'        , 'Entries'         , 'hdx'    )
hdy   = H1 ( 100, -  1,   1, 'hdy'  , 'dy (cm)'        , 'Entries'         , 'hdy'    )
hdz   = H1 ( 100, -  1,   1, 'hdz'  , 'dz (cm)'        , 'Entries'         , 'hdz'    )
hdr   = H1 ( 100,    0,   1, 'hdr'  , 'dr (cm)'        , 'Entries'         , 'hdr'    )
hdR   = H1 ( 100,    0,   1, 'hdR'  , 'dR (cm)'        , 'Entries'         , 'hdR'    )
hdEdR = H1 ( 100,    0, 0.3, 'hdEdR', 'dE/dR (MeV/cm)' , 'Entries'         , 'hdEdR'  )
hsx   = H1 ( 100, -  1,   1, 'hsx'  , 'sx (cm)'        , 'Entries'         , 'hsx'    )
hsy   = H1 ( 100, -  1,   1, 'hsy'  , 'sy (cm)'        , 'Entries'         , 'hsy'    )
hsz   = H1 ( 100, -  1,   1, 'hsz'  , 'sz (cm)'        , 'Entries'         , 'hsz'    )
hsr   = H1 ( 100, -  1,   1, 'hsr'  , 'sr (cm)'        , 'Entries'         , 'hsr'    )
hsR   = H1 ( 100, -  1,   1, 'hsr'  , 'sr (cm)'        , 'Entries'         , 'hsR'    )
hsrt  = H2P(  50,    0,  50,
             100,    0,   1, 'hsrt' , 'step'           , '#Deltar (cm)'    , 'hsrt'   )
hsxt  = H2P(  50,    0,  50,
             100, -  1,   1, 'hsxt' , 'step'           , '#Deltax (cm)'    , 'hsxt'   )
hsyt  = H2P(  50,    0,  50,
             100, -  1,   1, 'hsyt' , 'step'           , '#Deltay (cm)'    , 'hsyt'   )
htht  = H2P(  50,    0,  50,
             100,    0, pi5, 'htht' , 'step'           , '#theta'          , 'htht'   )
hthxt = H2P(  50,    0,  50,
             100, -pi5, pi5, 'hthxt', 'step'           , '#theta_{x}'      , 'hthxt'  )
hthyt = H2P(  50,    0,  50,
             100, -pi5, pi5, 'hthyt', 'step'           , '#theta_{y}'      , 'hthyt'  )


hxx   = H1 ( 100, - 1,   1, 'hxx'   , 'xx   correlation', 'Entries'        , 'hxx'    )
hxy   = H1 ( 100, - 1,   1, 'hxy'   , 'xy   correlation', 'Entries'        , 'hxy'    )
hxz   = H1 ( 100, - 1,   1, 'hxz'   , 'xz   correlation', 'Entries'        , 'hxz'    )
hxux  = H1 ( 100, - 1,   1, 'hxux'  , 'xux  correlation', 'Entries'        , 'hxux'   )
hxuy  = H1 ( 100, - 1,   1, 'hxuy'  , 'xuy  correlation', 'Entries'        , 'hxuy'   )
hxuz  = H1 ( 100, - 1,   1, 'hxuz'  , 'xuz  correlation', 'Entries'        , 'hxuz'   )
hxE   = H1 ( 100, - 1,   1, 'hxE'   , 'xE   correlation', 'Entries'        , 'hxE'    )
hxdE  = H1 ( 100, - 1,   1, 'hxdE'  , 'xdE  correlation', 'Entries'        , 'hxdE'   )
hxds  = H1 ( 100, - 1,   1, 'hxds'  , 'xds  correlation', 'Entries'        , 'hxds'   )

hyx   = H1 ( 100, - 1,   1, 'hyx'   , 'yx   correlation', 'Entries'         , 'hyx'   )
hyy   = H1 ( 100, - 1,   1, 'hyy'   , 'yy   correlation', 'Entries'         , 'hyy'   )
hyz   = H1 ( 100, - 1,   1, 'hyz'   , 'yz   correlation', 'Entries'         , 'hyz'   )
hyux  = H1 ( 100, - 1,   1, 'hyux'  , 'yux  correlation', 'Entries'         , 'hyux'  )
hyuy  = H1 ( 100, - 1,   1, 'hyuy'  , 'yuy  correlation', 'Entries'         , 'hyuy'  )
hyuz  = H1 ( 100, - 1,   1, 'hyuz'  , 'yuz  correlation', 'Entries'         , 'hyuz'  )
hyE   = H1 ( 100, - 1,   1, 'hyE'   , 'yE   correlation', 'Entries'         , 'hyE'   )
hydE  = H1 ( 100, - 1,   1, 'hydE'  , 'ydE  correlation', 'Entries'         , 'hydE'  )
hyds  = H1 ( 100, - 1,   1, 'hyds'  , 'yds  correlation', 'Entries'         , 'hyds'  )

hzx   = H1 ( 100, - 1,   1, 'hzx'   , 'zx   correlation', 'Entries'         , 'hzx'   )
hzy   = H1 ( 100, - 1,   1, 'hzy'   , 'zy   correlation', 'Entries'         , 'hzy'   )
hzz   = H1 ( 100, - 1,   1, 'hzz'   , 'zz   correlation', 'Entries'         , 'hzz'   )
hzux  = H1 ( 100, - 1,   1, 'hzux'  , 'zux  correlation', 'Entries'         , 'hzux'  )
hzuy  = H1 ( 100, - 1,   1, 'hzuy'  , 'zuy  correlation', 'Entries'         , 'hzuy'  )
hzuz  = H1 ( 100, - 1,   1, 'hzuz'  , 'zuz  correlation', 'Entries'         , 'hzuz'  )
hzE   = H1 ( 100, - 1,   1, 'hzE'   , 'zE   correlation', 'Entries'         , 'hzE'   )
hzdE  = H1 ( 100, - 1,   1, 'hzdE'  , 'zdE  correlation', 'Entries'         , 'hzdE'  )
hzds  = H1 ( 100, - 1,   1, 'hzds'  , 'zds  correlation', 'Entries'         , 'hzds'  )

huxx  = H1 ( 100, - 1,   1, 'huxx'  , 'uxx  correlation', 'Entries'         , 'huxx'  )
huxy  = H1 ( 100, - 1,   1, 'huxy'  , 'uxy  correlation', 'Entries'         , 'huxy'  )
huxz  = H1 ( 100, - 1,   1, 'huxz'  , 'uxz  correlation', 'Entries'         , 'huxz'  )
huxux = H1 ( 100, - 1,   1, 'huxux' , 'uxux correlation', 'Entries'         , 'huxux' )
huxuy = H1 ( 100, - 1,   1, 'huxuy' , 'uxuy correlation', 'Entries'         , 'huxuy' )
huxuz = H1 ( 100, - 1,   1, 'huxuz' , 'uxuz correlation', 'Entries'         , 'huxuz' )
huxE  = H1 ( 100, - 1,   1, 'huxE'  , 'uxE  correlation', 'Entries'         , 'huxE'  )
huxdE = H1 ( 100, - 1,   1, 'huxdE' , 'uxdE correlation', 'Entries'         , 'huxdE' )
huxds = H1 ( 100, - 1,   1, 'huxds' , 'uxds correlation', 'Entries'         , 'huxds' )

huyx  = H1 ( 100, - 1,   1, 'huyx'  , 'uyx  correlation', 'Entries'         , 'huyx'  )
huyy  = H1 ( 100, - 1,   1, 'huyy'  , 'uyy  correlation', 'Entries'         , 'huyy'  )
huyz  = H1 ( 100, - 1,   1, 'huyz'  , 'uyz  correlation', 'Entries'         , 'huyz'  )
huyux = H1 ( 100, - 1,   1, 'huyux' , 'uyux correlation', 'Entries'         , 'huyux' )
huyuy = H1 ( 100, - 1,   1, 'huyuy' , 'uyuy correlation', 'Entries'         , 'huyuy' )
huyuz = H1 ( 100, - 1,   1, 'huyuz' , 'uyuz correlation', 'Entries'         , 'huyuz' )
huyE  = H1 ( 100, - 1,   1, 'huyE'  , 'uyE  correlation', 'Entries'         , 'huyE'  )
huydE = H1 ( 100, - 1,   1, 'huydE' , 'uydE correlation', 'Entries'         , 'huydE' )
huyds = H1 ( 100, - 1,   1, 'huyds' , 'uyds correlation', 'Entries'         , 'huyds' )

huzx  = H1 ( 100, - 1,   1, 'huzx'  , 'uzx  correlation', 'Entries'         , 'huzx'  )
huzy  = H1 ( 100, - 1,   1, 'huzy'  , 'uzy  correlation', 'Entries'         , 'huzy'  )
huzz  = H1 ( 100, - 1,   1, 'huzz'  , 'uzz  correlation', 'Entries'         , 'huzz'  )
huzux = H1 ( 100, - 1,   1, 'huzux' , 'uzux correlation', 'Entries'         , 'huzux' )
huzuy = H1 ( 100, - 1,   1, 'huzuy' , 'uzuy correlation', 'Entries'         , 'huzuy' )
huzuz = H1 ( 100, - 1,   1, 'huzuz' , 'uzuz correlation', 'Entries'         , 'huzuz' )
huzE  = H1 ( 100, - 1,   1, 'huzE'  , 'uzE  correlation', 'Entries'         , 'huzE'  )
huzdE = H1 ( 100, - 1,   1, 'huzdE' , 'uzdE correlation', 'Entries'         , 'huzdE' )
huzds = H1 ( 100, - 1,   1, 'huzds' , 'uzds correlation', 'Entries'         , 'huzds' )

hEx   = H1 ( 100, - 1,   1, 'hEx'   , 'Ex   correlation', 'Entries'         , 'hEx'   )
hEy   = H1 ( 100, - 1,   1, 'hEy'   , 'Ey   correlation', 'Entries'         , 'hEy'   )
hEz   = H1 ( 100, - 1,   1, 'hEz'   , 'Ez   correlation', 'Entries'         , 'hEz'   )
hEux  = H1 ( 100, - 1,   1, 'hEux'  , 'Eux  correlation', 'Entries'         , 'hEux'  )
hEuy  = H1 ( 100, - 1,   1, 'hEuy'  , 'Euy  correlation', 'Entries'         , 'hEuy'  )
hEuz  = H1 ( 100, - 1,   1, 'hEuz'  , 'Euz  correlation', 'Entries'         , 'hEuz'  )
hEE   = H1 ( 100, - 1,   1, 'hEE'   , 'EE   correlation', 'Entries'         , 'hEE'   )
hEdE  = H1 ( 100, - 1,   1, 'hEdE'  , 'EdE  correlation', 'Entries'         , 'hEdE'  )
hEds  = H1 ( 100, - 1,   1, 'hEds'  , 'Eds  correlation', 'Entries'         , 'hEds'  )

hdEx  = H1 ( 100, - 1,   1, 'hdEx'  , 'dEx  correlation', 'Entries'         , 'hdEx'  )
hdEy  = H1 ( 100, - 1,   1, 'hdEy'  , 'dEy  correlation', 'Entries'         , 'hdEy'  )
hdEz  = H1 ( 100, - 1,   1, 'hdEz'  , 'dEz  correlation', 'Entries'         , 'hdEz'  )
hdEux = H1 ( 100, - 1,   1, 'hdEux' , 'dEux correlation', 'Entries'         , 'hdEux' )
hdEuy = H1 ( 100, - 1,   1, 'hdEuy' , 'dEuy correlation', 'Entries'         , 'hdEuy' )
hdEuz = H1 ( 100, - 1,   1, 'hdEuz' , 'dEuz correlation', 'Entries'         , 'hdEuz' )
hdEE  = H1 ( 100, - 1,   1, 'hdEE'  , 'dEE  correlation', 'Entries'         , 'hdEE'  )
hdEdE = H1 ( 100, - 1,   1, 'hdEdE' , 'dEdE correlation', 'Entries'         , 'hdEdE' )
hdEds = H1 ( 100, - 1,   1, 'hdEds' , 'dEds correlation', 'Entries'         , 'hdEds' )

hdsx   = H1 ( 100, - 1,   1, 'hdsx' , 'dsx  correlation', 'Entries'         , 'hdsx'  )
hdsy   = H1 ( 100, - 1,   1, 'hdsy' , 'dsy  correlation', 'Entries'         , 'hdsy'  )
hdsz   = H1 ( 100, - 1,   1, 'hdsz' , 'dsz  correlation', 'Entries'         , 'hdsz'  )
hdsux  = H1 ( 100, - 1,   1, 'hdsux', 'dsux correlation', 'Entries'         , 'hdsux' )
hdsuy  = H1 ( 100, - 1,   1, 'hdsuy', 'dsuy correlation', 'Entries'         , 'hdsuy' )
hdsuz  = H1 ( 100, - 1,   1, 'hdsuz', 'dsuz correlation', 'Entries'         , 'hdsuz' )
hdsE   = H1 ( 100, - 1,   1, 'hdsE' , 'dsE  correlation', 'Entries'         , 'hdsE'  )
hdsdE  = H1 ( 100, - 1,   1, 'hdsdE', 'dsdE correlation', 'Entries'         , 'hdsdE' )
hdsds  = H1 ( 100, - 1,   1, 'hdsds', 'dsds correlation', 'Entries'         , 'hdsds' )

hcorr  = HP2(   9,   0,   9,
                9,   0,   9, 'hcorr',                 '',        '',      '', 'hcorr' )
hcorr.GetXaxis().SetBinLabel( 1, 'x'  ); hcorr.GetYaxis().SetBinLabel( 1,  'x' )
hcorr.GetXaxis().SetBinLabel( 2, 'y'  ); hcorr.GetYaxis().SetBinLabel( 2,  'y' )
hcorr.GetXaxis().SetBinLabel( 3, 'z'  ); hcorr.GetYaxis().SetBinLabel( 3,  'z' )
hcorr.GetXaxis().SetBinLabel( 4, 'ux' ); hcorr.GetYaxis().SetBinLabel( 4, 'ux' )
hcorr.GetXaxis().SetBinLabel( 5, 'uy' ); hcorr.GetYaxis().SetBinLabel( 5, 'uy' )
hcorr.GetXaxis().SetBinLabel( 6, 'uz' ); hcorr.GetYaxis().SetBinLabel( 6, 'uz' )
hcorr.GetXaxis().SetBinLabel( 7,  'E' ); hcorr.GetYaxis().SetBinLabel( 7,  'E' )
hcorr.GetXaxis().SetBinLabel( 8, 'dE' ); hcorr.GetYaxis().SetBinLabel( 8, 'dE' )
hcorr.GetXaxis().SetBinLabel( 9, 'ds' ); hcorr.GetYaxis().SetBinLabel( 9, 'ds' )


print 'Number of files: ', Nfiles

print 'running...'
for i in xrange( Nfiles ):
    if not i % 1000:
        print '    file number',i
    
    data = x,y,z,zf,ux,uy,uz,E,dE,ds = ReadFile( filename.replace('*',str(i)), skip = 1 )
    
    hxx  .Fill( Correlation(  x,  x ) ); hcorr.Fill( 0, 0, Correlation(  x,  x ) )
    hxy  .Fill( Correlation(  x,  y ) ); hcorr.Fill( 0, 1, Correlation(  x,  y ) )
    hxz  .Fill( Correlation(  x,  z ) ); hcorr.Fill( 0, 2, Correlation(  x,  z ) )
    hxux .Fill( Correlation(  x, ux ) ); hcorr.Fill( 0, 3, Correlation(  x, ux ) )
    hxuy .Fill( Correlation(  x, uy ) ); hcorr.Fill( 0, 4, Correlation(  x, uy ) )
    hxuz .Fill( Correlation(  x, uz ) ); hcorr.Fill( 0, 5, Correlation(  x, uz ) )
    hxE  .Fill( Correlation(  x,  E ) ); hcorr.Fill( 0, 6, Correlation(  x,  E ) )
    hxdE .Fill( Correlation(  x, dE ) ); hcorr.Fill( 0, 7, Correlation(  x, dE ) )
    hxds .Fill( Correlation(  x, ds ) ); hcorr.Fill( 0, 8, Correlation(  x, ds ) )
    
    hyx  .Fill( Correlation(  y,  x ) ); hcorr.Fill( 1, 0, Correlation(  y,  x ) )
    hyy  .Fill( Correlation(  y,  y ) ); hcorr.Fill( 1, 1, Correlation(  y,  y ) )
    hyz  .Fill( Correlation(  y,  z ) ); hcorr.Fill( 1, 2, Correlation(  y,  z ) )
    hyux .Fill( Correlation(  y, ux ) ); hcorr.Fill( 1, 3, Correlation(  y, ux ) )
    hyuy .Fill( Correlation(  y, uy ) ); hcorr.Fill( 1, 4, Correlation(  y, uy ) )
    hyuz .Fill( Correlation(  y, uz ) ); hcorr.Fill( 1, 5, Correlation(  y, uz ) )
    hyE  .Fill( Correlation(  y,  E ) ); hcorr.Fill( 1, 6, Correlation(  y,  E ) )
    hydE .Fill( Correlation(  y, dE ) ); hcorr.Fill( 1, 7, Correlation(  y, dE ) )
    hyds .Fill( Correlation(  y, ds ) ); hcorr.Fill( 1, 8, Correlation(  y, ds ) )
    
    hzx  .Fill( Correlation(  z,  x ) ); hcorr.Fill( 2, 0, Correlation(  z,  x ) )
    hzy  .Fill( Correlation(  z,  y ) ); hcorr.Fill( 2, 1, Correlation(  z,  y ) )
    hzz  .Fill( Correlation(  z,  z ) ); hcorr.Fill( 2, 2, Correlation(  z,  z ) )
    hzux .Fill( Correlation(  z, ux ) ); hcorr.Fill( 2, 3, Correlation(  z, ux ) )
    hzuy .Fill( Correlation(  z, uy ) ); hcorr.Fill( 2, 4, Correlation(  z, uy ) )
    hzuz .Fill( Correlation(  z, uz ) ); hcorr.Fill( 2, 5, Correlation(  z, uz ) )
    hzE  .Fill( Correlation(  z,  E ) ); hcorr.Fill( 2, 6, Correlation(  z,  E ) )
    hzdE .Fill( Correlation(  z, dE ) ); hcorr.Fill( 2, 7, Correlation(  z, dE ) )
    hzds .Fill( Correlation(  z, ds ) ); hcorr.Fill( 2, 8, Correlation(  z, ds ) )
    
    huxx .Fill( Correlation( ux,  x ) ); hcorr.Fill( 3, 0, Correlation( ux,  x ) )
    huxy .Fill( Correlation( ux,  y ) ); hcorr.Fill( 3, 1, Correlation( ux,  y ) )
    huxz .Fill( Correlation( ux,  z ) ); hcorr.Fill( 3, 2, Correlation( ux,  z ) )
    huxux.Fill( Correlation( ux, ux ) ); hcorr.Fill( 3, 3, Correlation( ux, ux ) )
    huxuy.Fill( Correlation( ux, uy ) ); hcorr.Fill( 3, 4, Correlation( ux, uy ) )
    huxuz.Fill( Correlation( ux, uz ) ); hcorr.Fill( 3, 5, Correlation( ux, uz ) )
    huxE .Fill( Correlation( ux,  E ) ); hcorr.Fill( 3, 6, Correlation( ux,  E ) )
    huxdE.Fill( Correlation( ux, dE ) ); hcorr.Fill( 3, 7, Correlation( ux, dE ) )
    huxds.Fill( Correlation( ux, ds ) ); hcorr.Fill( 3, 8, Correlation( ux, ds ) )
    
    huyx .Fill( Correlation( uy,  x ) ); hcorr.Fill( 4, 0, Correlation( uy,  x ) )
    huyy .Fill( Correlation( uy,  y ) ); hcorr.Fill( 4, 1, Correlation( uy,  y ) )
    huyz .Fill( Correlation( uy,  z ) ); hcorr.Fill( 4, 2, Correlation( uy,  z ) )
    huyux.Fill( Correlation( uy, ux ) ); hcorr.Fill( 4, 3, Correlation( uy, ux ) )
    huyuy.Fill( Correlation( uy, uy ) ); hcorr.Fill( 4, 4, Correlation( uy, uy ) )
    huyuz.Fill( Correlation( uy, uz ) ); hcorr.Fill( 4, 5, Correlation( uy, uz ) )
    huyE .Fill( Correlation( uy,  E ) ); hcorr.Fill( 4, 6, Correlation( uy,  E ) )
    huydE.Fill( Correlation( uy, dE ) ); hcorr.Fill( 4, 7, Correlation( uy, dE ) )
    huyds.Fill( Correlation( uy, ds ) ); hcorr.Fill( 4, 8, Correlation( uy, ds ) )
    
    huzx .Fill( Correlation( uz,  x ) ); hcorr.Fill( 5, 0, Correlation( uz,  x ) )
    huzy .Fill( Correlation( uz,  y ) ); hcorr.Fill( 5, 1, Correlation( uz,  y ) )
    huzz .Fill( Correlation( uz,  z ) ); hcorr.Fill( 5, 2, Correlation( uz,  z ) )
    huzux.Fill( Correlation( uz, ux ) ); hcorr.Fill( 5, 3, Correlation( uz, ux ) )
    huzuy.Fill( Correlation( uz, uy ) ); hcorr.Fill( 5, 4, Correlation( uz, uy ) )
    huzuz.Fill( Correlation( uz, uz ) ); hcorr.Fill( 5, 5, Correlation( uz, uz ) )
    huzE .Fill( Correlation( uz,  E ) ); hcorr.Fill( 5, 6, Correlation( uz,  E ) )
    huzdE.Fill( Correlation( uz, dE ) ); hcorr.Fill( 5, 7, Correlation( uz, dE ) )
    huzds.Fill( Correlation( uz, ds ) ); hcorr.Fill( 5, 8, Correlation( uz, ds ) )
    
    hEx  .Fill( Correlation(  E,  x ) ); hcorr.Fill( 6, 0, Correlation(  E,  x ) )
    hEy  .Fill( Correlation(  E,  y ) ); hcorr.Fill( 6, 1, Correlation(  E,  y ) )
    hEz  .Fill( Correlation(  E,  z ) ); hcorr.Fill( 6, 2, Correlation(  E,  z ) )
    hEux .Fill( Correlation(  E, ux ) ); hcorr.Fill( 6, 3, Correlation(  E, ux ) )
    hEuy .Fill( Correlation(  E, uy ) ); hcorr.Fill( 6, 4, Correlation(  E, uy ) )
    hEuz .Fill( Correlation(  E, uz ) ); hcorr.Fill( 6, 5, Correlation(  E, uz ) )
    hEE  .Fill( Correlation(  E,  E ) ); hcorr.Fill( 6, 6, Correlation(  E,  E ) )
    hEdE .Fill( Correlation(  E, dE ) ); hcorr.Fill( 6, 7, Correlation(  E, dE ) )
    hEds .Fill( Correlation(  E, ds ) ); hcorr.Fill( 6, 8, Correlation(  E, ds ) )
    
    hdEx .Fill( Correlation( dE,  x ) ); hcorr.Fill( 7, 0, Correlation( dE,  x ) )
    hdEy .Fill( Correlation( dE,  y ) ); hcorr.Fill( 7, 1, Correlation( dE,  y ) )
    hdEz .Fill( Correlation( dE,  z ) ); hcorr.Fill( 7, 2, Correlation( dE,  z ) )
    hdEux.Fill( Correlation( dE, ux ) ); hcorr.Fill( 7, 3, Correlation( dE, ux ) )
    hdEuy.Fill( Correlation( dE, uy ) ); hcorr.Fill( 7, 4, Correlation( dE, uy ) )
    hdEuz.Fill( Correlation( dE, uz ) ); hcorr.Fill( 7, 5, Correlation( dE, uz ) )
    hdEE .Fill( Correlation( dE,  E ) ); hcorr.Fill( 7, 6, Correlation( dE,  E ) )
    hdEdE.Fill( Correlation( dE, dE ) ); hcorr.Fill( 7, 7, Correlation( dE, dE ) )
    hdEds.Fill( Correlation( dE, ds ) ); hcorr.Fill( 7, 8, Correlation( dE, ds ) )

    hdsx .Fill( Correlation( ds,  x ) ); hcorr.Fill( 8, 0, Correlation( ds,  x ) )
    hdsy .Fill( Correlation( ds,  y ) ); hcorr.Fill( 8, 1, Correlation( ds,  y ) )
    hdsz .Fill( Correlation( ds,  z ) ); hcorr.Fill( 8, 2, Correlation( ds,  z ) )
    hdsux.Fill( Correlation( ds, ux ) ); hcorr.Fill( 8, 3, Correlation( ds, ux ) )
    hdsuy.Fill( Correlation( ds, uy ) ); hcorr.Fill( 8, 4, Correlation( ds, uy ) )
    hdsuz.Fill( Correlation( ds, uz ) ); hcorr.Fill( 8, 5, Correlation( ds, uz ) )
    hdsE .Fill( Correlation( ds,  E ) ); hcorr.Fill( 8, 6, Correlation( ds,  E ) )
    hdsdE.Fill( Correlation( ds, dE ) ); hcorr.Fill( 8, 7, Correlation( ds, dE ) )
    hdsds.Fill( Correlation( ds, ds ) ); hcorr.Fill( 8, 8, Correlation( ds, ds ) )


    data = zip( *data )
    nsteps = len( data )
    
    for j in range( 1, nsteps - 1 ):
        x,y,z,zf,ux,uy,uz,E,dE,ds = data[j]
        dx = data[j+1][0] - x
        dy = data[j+1][1] - y
        dz = data[j+1][2] - z
        dr = ( dx**2 + dy**2 )**.5
        dR = ( dr**2 + dz**2 )**.5
        xp = x + ux * ds
        yp = y + uy * ds
        zp = z + uz * ds
        Dx = xp - x - dx
        Dy = yp - y - dy
        Dz = zp - z - dz
        Dr = ( Dx**2 + Dy**2 )**.5
        DR = ( Dr**2 + Dz**2 )**.5
        xt = atan(Dx/ds)
        yt = atan(Dy/ds)
        rt = atan(DR/ds)

        hx   .Fill(  x )
        hy   .Fill(  y )
        hz   .Fill(  z )
        hux  .Fill( ux )
        huy  .Fill( uy )
        huz  .Fill( uz )
        hE   .Fill(  E )
        hdE  .Fill( dE )
        hds  .Fill( ds )
        hdx  .Fill( dx )
        hdy  .Fill( dy )
        hdz  .Fill( dz )
        hdr  .Fill( dr )
        hdR  .Fill( dR )
        hdEdR.Fill( dE / dR )
        hsx  .Fill( Dx )
        hsy  .Fill( Dy )
        hsz  .Fill( Dz )
        hsr  .Fill( Dr )
        hsrt .Fill(  j , Dr )
        hsxt .Fill(  j , Dx )
        hsyt .Fill(  j , Dy )
        htht .Fill(  j , rt )
        hthxt.Fill(  j , xt )
        hthyt.Fill(  j , yt )

a = PutInCanvas([hx,hy,hz,hux,huy,huz,hE,hdE,hds],None,3,2)
b = PutInCanvas([hdx,hdy,hdz,hdr,hdR,hdEdR], None, 3, 2)
c = PutInCanvas([hdsx,hdsy,hdsz,hdsux,hdsuy,hdsuz,hdsE,hdsdE,hdsds,
                 hdEx,hdEy,hdEz,hdEux,hdEuy,hdEuz,hdEE,hdEdE,hdEds,
                 hEx,hEy,hEz,hEux,hEuy,hEuz,hEE,hEdE,hEds,
                 huzx,huzy,huzz,huzux,huzuy,huzuz,huzE,huzdE,huzds,
                 huyx,huyy,huyz,huyux,huyuy,huyuz,huyE,huydE,huyds,
                 huxx,huxy,huxz,huxux,huxuy,huxuz,huxE,huxdE,huxds,
                 hzx,hzy,hzz,hzux,hzuy,hzuz,hzE,hzdE,hzds,
                 hyx,hyy,hyz,hyux,hyuy,hyuz,hyE,hydE,hyds,
                 hxx,hxy,hxz,hxux,hxuy,hxuz,hxE,hxdE,hxds],None,9,9)
d = PutInCanvas([hcorr],['textzcol'])
e = PutInCanvas([hsx,hsy,hsz,hsr])
f = PutInCanvas([hsxt,hsyt,hsrt,htht,hthxt,hthyt])