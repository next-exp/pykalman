'''
    Algorithms that read and ellaborate Josh's tracks.
'''

from ialex import IAlg, IReader
#from KalmanFMWK import KalmanMeasurement, KalmanNode, KalmanTrack

class JoshReader( IReader ):
    
    def open( self ):
        self.file = open( self.name, 'r' )

    def eof( self ):
        return self.line is ''

    def read( self ):
        return [ map( float, line.split() ) for line in self.file.readlines() ]

    def close( self ):
        self.file.close()
#
#
#class JoshTrackBuilder( IAlg ):
#    def __init__( self, name = 'JoshTrackBuilder' ):
#        IAlg.__init__( self, name )
#    
#    def define( self ):
#        self.inputpath = 'JoshReader'
#
#    def execute( self ):
#        data = self.evt.get(self.inputpath)
#        track = KalmanTrack()
#        for seg_id, k,x0,y0,z0,p1p,p2p,p3p,p4p,chi2p,p1f,p2f,p3f,p4f,chi2f,cfxy,cftxy,edep in data:
#            node = KalmanNode( step         = k
#                               hit          = KalmanMeasurement( Array.Vector( x0, y0, z0 ) )
#                               pred_state   = KalmanMeasurement( Array.Vector( p1p, p2p, p3p, p4p ) ),
#                               filt_state   = KalmanMeasurement( Array.Vector( p1f, p2f, p3f, p4f ) ),
#                               smooth_state = KalmanMeasurement()
#                               pred_resid   = KalmanMeasurement(),
#                               filt_resid   = KalmanMeasurement(),
#                               smooth_resid = KalmanMeasurement(),
#                               chi2         = chi2f, cumchi2 = -1 )
#            track.AddNode( node )
#        
#        return True
#
#    def finalize( self ):
#        return
