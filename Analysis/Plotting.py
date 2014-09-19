'''
    Plotting.py

    Description:
    Some functions and definitions in order to make the analysis more comfortable.
    
    Author:  Gonzalo
    Date:   20/07/2014
    
    Last update 20/07/2014 <- Gonzalo
'''

from array import array
from numpy import array as narray
from ROOT import *

class H2D(list):
    '''Contains methods to apply histograms methods both to scatter and profile histograms.'''
    
    def Fill(self,x,y):
        map(lambda h: h.Fill(x,y), self)
    
    def Draw(self,who='both',opt=''):
        if who=='2':
            self[0].Draw(opt)
            return None
        elif who=='P':
            self[1].Draw(opt)
            return None
        else:
            self[1].SetLineColor(2)
            map(lambda x: x.Draw('SAME'), self)
        return None
    
    def Write(self):
        map(lambda x: x.Write(), self)

class H3D(H2D):
    '''Contains methods to apply histograms methods both to scatter and profile histograms.'''
    
    def Fill(self,X,Y,Z):
        map(lambda x: x.Fill(X,Y,Z), self)
    
    def Draw(self,who='P',opt=''):
        if who=='3':
            self[0].Draw(opt)
            return None
        elif who=='P':
            self[1].Draw(opt)
            return None
        else:
            raise ValueError('Wrong argument')
        return None


def H1(nbins,min,max,title='',xlabel='',ylabel='',name=None):
    '''1D histogram.'''
    
    if not name:
        name = title
    
    title += ';' + xlabel + ';' + ylabel
    
    return TH1F(name,title,nbins,min,max)

def H2(nbinsX,minX,maxX,nbinsY,minY,maxY,title='',xlabel='',ylabel='',name=None):
    '''2D scatter histogram.'''
    
    if not name:
        name = title
    
    title += ';' + xlabel + ';' + ylabel
    
    return TH2F(name,title,nbinsX,minX,maxX,nbinsY,minY,maxY)

def HP(nbins,min,max,title='',xlabel='',ylabel='',name=None):
    '''1D profile histogram.'''
    
    if not name:
        name = title
    
    title += ';' + xlabel + ';' + ylabel
    
    return TProfile(name,title,nbins,min,max)

def H2P(nbinsX,minX,maxX,nbinsY,minY,maxY,title='',xlabel='',ylabel='',name=''):
    '''Both 2D scatter and its profile histograms into the H2D class.'''
    
    return H2D([H2(nbinsX,minX,maxX,nbinsY,minY,maxY,title,xlabel,ylabel,name           ),
                HP(nbinsX,minX,maxX,                 title,xlabel,ylabel,name+'_profile') ])

def H3(nbinsX,minX,maxX,nbinsY,minY,maxY,nbinsZ,minZ,maxZ,title='',xlabel='',ylabel='',zlabel='',name=None):
    '''3D scatter histogram.'''
    
    if not name:
        name = title
    
    title += ';' + xlabel + ';' + ylabel + ';' + zlabel
    
    return TH3F(name,title,nbinsX,minX,maxX,nbinsY,minY,maxY,nbinsZ,minZ,maxZ)

def HP2(nbinsX,minX,maxX,nbinsY,minY,maxY,title='',xlabel='',ylabel='',zlabel='',name=None):
    '''2D profile histogram.'''
    
    if not name:
        name = title
    
    title += ';' + xlabel + ';' + ylabel + ';' + zlabel
    
    return TProfile2D(name,title,nbinsX,minX,maxX,nbinsY,minY,maxY)

def H3P(nbinsX,minX,maxX,nbinsY,minY,maxY,nbinsZ,minZ,maxZ,title='',xlabel='',ylabel='',zlabel='',name=''):
    '''Both 3D scatter and its profile histograms into the H3D class.'''
    
    return H3D([ H3(nbinsX,minX,maxX,nbinsY,minY,maxY,nbinsZ,minZ,maxZ,title,xlabel,ylabel,zlabel,name ),
                HP2(nbinsX,minX,maxX,nbinsY,minY,maxY,                 title,xlabel,ylabel,zlabel,name+'_profile' ) ])

Colors = [  4,   2,   1,  53,   8,
           95,  14,   7, 205,   3,
          218,   6,  46,  40,   5] * 100

canvasorganization     = {}
canvasorganization[1]  = ( 1, 1 )
canvasorganization[2]  = ( 2, 1 )
canvasorganization[3]  = ( 3, 1 )
canvasorganization[4]  = ( 2, 2 )
canvasorganization[5]  = ( 2, 3 )
canvasorganization[6]  = ( 2, 3 )
canvasorganization[7]  = ( 2, 4 )
canvasorganization[8]  = ( 2, 4 )
canvasorganization[9]  = ( 3, 3 )
canvasorganization[10] = ( 3, 4 )
canvasorganization[11] = ( 3, 4 )
canvasorganization[12] = ( 3, 4 )
canvasorganization[13] = ( 4, 4 )
canvasorganization[14] = ( 4, 4 )
canvasorganization[15] = ( 4, 4 )
canvasorganization[16] = ( 4, 4 )
canvasorganization[17] = ( 4, 5 )
canvasorganization[18] = ( 4, 5 )
canvasorganization[19] = ( 4, 5 )
canvasorganization[20] = ( 4, 5 )
canvasorganization[21] = ( 5, 5 )
canvasorganization[22] = ( 5, 5 )
canvasorganization[23] = ( 5, 5 )
canvasorganization[24] = ( 5, 5 )

def printcolortable():
    ''' Prints the color table that which is used when plots are merged.'''
    
    COLORS = Colors[:15]
    NofC = len(COLORS)
    
    c = TCanvas()
    h = [ TH1I( 'colortable' + str(i), 'Color table', NofC, 0, NofC ) for i in range(NofC) ]
    
    gStyle.SetOptStat('')
    map( lambda i, c: h[i].SetFillColor(c)                     , range(NofC), COLORS )
    map( lambda hist: hist.GetXaxis().SetNdivisions(NofC,False), h )
    map( lambda hist: hist.GetYaxis().SetLabelColor(0), h )
    map( lambda hist: hist.GetYaxis().SetTickLength(0), h )
    map( lambda hist: hist.SetMaximum(1), h )
    
    for i in range(NofC):
        h[i].Fill(i)
        h[i].Draw('same')
    
    return c,h

def Graph( x, y, xaxis='', yaxis='',title='', markerstyle = None, markersize = None ):
    ''' Returns a TGraph made of x,y data with optional titles.'''
    
    if not len(x) == len(y):
        raise ValueError('Both data lists must have the same length.')
    
    graph = TGraph( len(x), array('f',x), array('f',y) )
    
    if xaxis:
        graph.GetXaxis().SetTitle( xaxis )
    if yaxis:
        graph.GetYaxis().SetTitle( yaxis )
    if title:
        graph.SetTitle( title )
    if markerstyle:
        graph.SetMarkerStyle( markerstyle )
    if markersize:
        graph.SetMarkerSize( markersize )
    
    return graph

def ErrorsGraph( x, y, ex, ey, xaxis='', yaxis='',title='' ):
    ''' Returns a TGraphErrors made of x,y data with optional titles.'''
    
    if ex is None:
        ex = [0.] * len(x)
    if ey is None:
        ey = [0.] * len(y)

    if not len(x) == len(y) == len(ex) == len(ey):
        raise ValueError('Both data lists must have the same length.')
    
    graph = TGraphErrors( len(x), array('f',x), array('f',y), array('f',ex), array('f',ey) )
    
    if xaxis:
        graph.GetXaxis().SetTitle( xaxis )
    if yaxis:
        graph.GetYaxis().SetTitle( yaxis )
    if title:
        graph.SetTitle( title )
    
    return graph

def GetMax( histogram ):
    ''' Returns the value of the bin with the greatest value.'''
    
    return histogram.GetBinContent( histogram.GetMaximumBin() )

def Norm2integral( histogram ):
    ''' Scales the histogram so that the integral is 1.'''
    
    histogram.Scale( 1./ histogram.Integral() )
    return histogram

def Norm2maximum( histogram ):
    ''' Scales the histogram so that the greatest bin is setted to 1.'''
    
    histogram.Scale( 1./ GetMax( histogram ) )
    return histogram

def Merge( histograms, xaxis = '', yaxis = '', title = '' ):
    ''' Draws in the same canvas all the histograms passed. Returns the canvas.'''
    
    c = TCanvas()
    
    for i in range( len(histograms) ):
        h = histograms[i]
        h.SetLineColor( Colors[i] )
        h.SetMarkerColor( Colors[i] )
        h.Draw( 'SAME' )
    
    histograms[0].SetMaximum( max( map( GetMax, histograms ) ) )
    
    if xaxis:
        histograms[0].GetXaxis().SetTitle( xaxis )
    if yaxis:
        histograms[0].GetYaxis().SetTitle( yaxis )
    if title:
        histograms[0].SetTitle( title )
    
    return c

def Gethisto( file, hname ):
    ''' Gets an histogram from a file. It accepts the name of the file or the file itself. Use this last option when getting a lot of histograms from the same file.'''
    
    if isinstance( file, str ):
        file = TFile( file )
    
    h = file.Get( hname )
    if isinstance( h, TProfile ):
        h.SetErrorOption('')
    
    return h

def Addprofile( h2, hp, color = 2, width = 3 ):
    ''' Merges the scatter plot with its profile. Returns the canvas.'''
    
    c = TCanvas()
    h2.Draw()
    hp.Draw('SAME')
    hp.SetLineColor( color )
    hp.SetLineWidth( width )
    
    return c

def Sumhistos( *hlist ):
    ''' Sums the histograms.'''
    
    if not isinstance( hlist, (list,tuple) ):
        hlist = [ hlist ]
    
    H = hlist[0]
    for h in hlist[1:]:
        H.Add( h )

def GoodLooking( histogram, color = 1, width = 2, fill = None ):
    ''' Sets the usual attributer to the histogram for a fancy presentation.'''
    
    histogram.SetLineColor( color )
    histogram.SetLineWidth( width )
    if fill:
        histogram.SetFillColor( fill )

def MakeH1( data, title = 'histo', nbins = 100 ):
    ''' Returns the distribution of data.'''
    
    Data = sorted( data )
    if isinstance( data[0], (tuple) ):
        nbins = len( Data )
        MIN   = Data[0][0]
        MAX   = Data[-1][0]
        MAX  += ( MAX - MIN ) / ( nbins - 1 )
        histo = TH1F( title, title, nbins, MIN, MAX )
        [ histo.SetBinContent( i + 1, data[i][1] ) for i in range( nbins ) ]
    else:
        MIN  = Data[ 0]
        MAX  = Data[-1]
        MAX += ( MAX - MIN ) / ( nbins - 1 )
        histo = TH1F( title, title, nbins, MIN, MAX )
        map( histo.Fill, data )
    
    return histo

def PutInCanvas( objects, Options = None, nhorizontal = None, nvertical = None ):
    if nhorizontal is None and nvertical is None:
        nhorizontal, nvertical = canvasorganization[ len(objects) ]
    
    if Options is None:
        Options = [''] * len( objects )
    
    c = TCanvas()
    c.Divide( nhorizontal, nvertical )
    
    for i in range(len(objects)):
        c.cd(i+1)
        objects[i].Draw( Options[i] )
    
    return c

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
