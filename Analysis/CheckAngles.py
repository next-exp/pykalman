'''
    This file shows how the dispersion in the scattering angle affects to the spherical angles Theta and Phi.
'''


import ROOT
import math
import Plots
pi = math.pi
ROOT.gStyle.SetOptStat('')

R = ROOT.TRandom3(0)
U = R.Uniform
G = R.Gaus
V = lambda x: (2*U()-1)*x

z      = 1.0
sigmas = 0.2 #small
sigmal = 0.8 #large

hthetas = ROOT.TH1F( 'smalltheta', '#theta_{x} = #theta_{y};#theta_{x} = #theta_{y};#Entries', 250, -pi,   pi )
hthetal = ROOT.TH1F( 'largetheta', '#theta_{x} = #theta_{y};#theta_{x} = #theta_{y};#Entries', 250, -pi,   pi )
htans   = ROOT.TH1F( 'tansmall'  , 'tan x = tan y;tan x = tan y;Entries'                     , 250, -20,   20 )
htanl   = ROOT.TH1F( 'tanlarge'  , 'tan x = tan y;tan x = tan y;Entries'                     , 250, -20,   20 )
hrhos   = ROOT.TH1F( 'rhosmall'  , '#rho;#rho (a.u.);Entries'                                , 250,   0,   20 )
hrhol   = ROOT.TH1F( 'rholarge'  , '#rho;#rho (a.u.);Entries'                                , 250,   0,   20 )
hThetas = ROOT.TH1F( 'Thetasmall','#Theta;#Theta;Entries'                                    , 250,  0,  pi/2 )
hThetal = ROOT.TH1F( 'Thetalarge','#Theta;#Theta;Entries'                                    , 250,  0,  pi/2 )
hphis   = ROOT.TH1F( 'phismall'  ,'#Phi;#Phi;Entries'                                        , 250,  0,  2*pi )
hphil   = ROOT.TH1F( 'philarge'  ,'#Phi;#Phi;Entries'                                        , 250,  0,  2*pi )
hxys    = ROOT.TH2F( 'xysmall'   , 'xy;x (a.u.);y (a.u.)'                                    , 250, -20,   20 , 250, -20, 20 )
hxyl    = ROOT.TH2F( 'xylarge'   , 'xy;x (a.u.);y (a.u.)'                                    , 250, -20,   20 , 250, -20, 20 )


hthetas.SetMinimum(0)
hthetal.SetMinimum(0)
htans  .SetMinimum(0)
htanl  .SetMinimum(0)
hrhos  .SetMinimum(0)
hrhol  .SetMinimum(0)
hxys   .SetMinimum(0)
hxyl   .SetMinimum(0)
hThetas.SetMinimum(0)
hThetal.SetMinimum(0)
hphis  .SetMinimum(0)
hphil  .SetMinimum(0)

for i in xrange(200000):

    if i % 2:
        thetax = G( 0, sigmas )
        thetay = G( 0, sigmas )
    else:
        thetax = G( 0, sigmal )
        thetay = G( 0, sigmal )

    tanx = math.tan( thetax )
    tany = math.tan( thetay )

    x = z * tanx
    y = z * tany
    z = z

    rho = (   x**2 + y**2 )**.5
    r   = ( rho**2 + z**2 )**.5

    theta = math.acos( z / r )
    phi   = math.atan( y / x )
    phi  += pi if x<0 else 2*pi if y<0 else 0

    if i % 2:
        hthetas.Fill( thetax )
        hthetas.Fill( thetay )
        htans  .Fill(   tanx )
        htans  .Fill(   tany )
        hrhos  .Fill(   rho  )
        hThetas.Fill( theta  )
        hphis  .Fill(   phi  )
        hxys   .Fill(   x, y )

    else:
        hthetal.Fill( thetax )
        hthetal.Fill( thetay )
        htanl  .Fill(   tanx )
        htanl  .Fill(   tany )
        hrhol  .Fill(   rho  )
        hThetal.Fill( theta  )
        hphil  .Fill(   phi  )
        hxyl   .Fill(   x, y )


canvas1 = Plots.PutInCanvas( [ hthetas, hthetal,
                               htans  , htanl  ,
                               hrhos  , hrhol  ,
                               hThetas, hThetal,
                               hphis  , hphil  ,
                               hxys   , hxyl   ], ['']*10 + ['zcol'] * 2, 2, 6 )


#canvas1 = PutInCanvas([hthetax,hthetay,htanx,htany,htheta,hphi,hrho,hxy],['']*7 + ['zcol'],2,4)