from __future__ import division
from os import listdir

def rangelen( list ):
    return range(len(list))

def Mean( x, w = None ):
    '''
        Return the mean of x with weigths w.
    '''
    return sum([ xi*wi for xi,wi in zip(x,w) ])/sum(w) if w else sum(x)/len(x)

def Covariance( x, y, w = None ):
    '''
        Return the covariance between data in x and y.
    '''
    xmean = Mean( x, w )
    ymean = Mean( y, w )
    if w:
        sumw  = sum(w)
        sumw2 = sum([ wi**2 for wi in w ])
        norm  = sumw / ( sumw**2 - sumw2 )
        return norm * sum( [ wi * (xi-xmean) * (yi-ymean) for xi,yi,wi in zip(x,y,w) ] )
    else:
        return sum( [ (xi-xmean) * (yi-ymean) for xi,yi in zip(x,y) ] ) / len(x)

def RMS( x, w = None ):
    return Covariance( x, x, w )**.5

def Correlation( x, y, w = None ):
    return Covariance( x, y, w ) / ( RMS( x, w ) * RMS( y, w ) )

def ReadFile( filename, separator = ' ', skip = 0 ):
    ''' Reads data from a file and returns a list of lists of data. If there is more than one column of data, you can specify the separator between them which is a space by default. Also, in order to convert this data from strings to numbers you can specify the type of numbers, which is taken like floats by default, but you can choose also integers.'''
    
    datain = [ map( float, line.split(separator) ) for i,line in enumerate( open( filename, 'r') ) if i>=skip ]

    return [ [ datain[j][i] for j in rangelen(datain) ] for i in rangelen( datain[0] ) ]

def GetNFiles( name ):
    return len(listdir( name )) - 1
    