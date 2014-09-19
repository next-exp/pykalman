'''
    Kalman stuff
'''
import copy
import Array
#
#class KalmanHit:
#    
#    def __init__( self, x, y, z, E, CovarianceMatrix ):
#        Array.Vector4.__init__( self, E, x, y, z )
#        self.CovarianceMatrix = CovarianceMatrix
#
#    def __str__( self ):
#        return '''Hit = {0}
#                  Cov = 
#                  {1}
#                  '''.format( Array.Vector4.__str__(self), self.CovarianceMatrix)
#

class KalmanMeasurement:
    '''
        Represents a state of the Kalman Filter.
    '''

    def __init__( self, state = Array.Vector(), CovarianceMatrix = Array.Matrix() ):
        self.state = state
        self.CovarianceMatrix = CovarianceMatrix

    def __str__( self ):
        return '''state = {0}
                  Cov =
                  {1}
                  '''.format( self.state, self.CovarianceMatrix )

class KalmanNode:
    '''
        Contains all the information about a given step.
    '''
    
    def __init__( self, step = 0, hit = KalmanMeasurement(),
                  pred_state = KalmanMeasurement(), filt_state = KalmanMeasurement(), smooth_state = KalmanMeasurement(),
                  pred_resid = KalmanMeasurement(), filt_resid = KalmanMeasurement(), smooth_resid = KalmanMeasurement(),
                  chi2       = -1, cumchi2 = -1 ):
        self.step = step
        self.hit = hit
        self.pred_state   = pred_state
        self.filt_state   = filt_state
        self.smooth_state = smooth_state
        self.pred_resid   = pred_resid
        self.filt_resid   = filt_resid
        self.smooth_resid = smooth_resid
        self.chi2         = chi2
        self.cumchi2      = cumchi2

    def __str__( self ):
        return '''step number {0}
                  hit:
                  {1}
                  predicted state:
                  {2}
                  filtered state:
                  {3}
                  smoothed state:
                  {4}
                  predicted residuals:
                  {5}
                  filtered residuals:
                  {6}
                  smoothed residuals:
                  {7}
                  chi2 of this node:
                  {8}
                  cumulative chi2 at this node:
                  {9}
                  '''.format( self.step, self.hit, self.pred_state, self.filt_state, self.smooth_state,
                                                   self.pred_resid, self.filt_resid, self.smooth_resid,
                              self.chi2, self.cumchi2 )

class KalmanTrack:
    '''
        A set of KalmanNodes.
    '''
    def __init__( self, nodes = list() ):
        self.nodes = copy.deepcopy(nodes)
        self.nnodes = len( self.nodes )

    def AddNode( self, node ):
        self.nodes.append( copy.copy( node ) )
        self.nnodes += 1
    
    def GetNode( self, index ):
        return self.nodes[index]

    def __str__( self ):
        return 'Number of steps: {0} \n\n'.format( self.nnodes ) + '\n\n'.join( map( str, self.nodes ) )

class KalmanFilter:
    '''
        Abstract implementation of the Kalman Filter.
    '''
    def __init__( self, name = 'KalmanFilter' ):
        self.name             = name
        self.Track            = KalmanTrack()
        self.HaveInitialState = False
        self.HaveMeasurements = False
        self.guess            = None

    def SetInitialState( self, state = KalmanMeasurement() ):
        self.Track.GetNode(0).pred_state = state
        self.Track.GetNode(0).filt_state = state
        self.ndim  = len(state)
        self.HaveInitialState = True
    
    def SetInitialGuess( self, guess = None ):
        self.guess = guess if guess else self._ComputeInitialGuess()
    
    def SetMeasurements( self, hits ):
        for i,hit in enumerate(hits):
            self.Track.AddNode( KalmanNode(step = i, hit = hit) )
        self.HaveMeasurements = True

    def _ComputeInitialGuess( self ):
        #To be implemented
        pass

    def TransportMatrix( self, index ):
        '''
            To be coded in the particular implementation.
        '''
        return
    
    def MeasurementMatrix( self, index ):
        '''
            To be coded in the particular implementation.
        '''
        return
    
    def NoiseMatrix( self, index ):
        '''
            To be coded in the particular implementation.
        '''
        return

    def _NewNode( self, index ):
        
        self.prev_node  = self.Track.GetNode( index - 1 )
        self.this_node  = self.Track.GetNode( index )
        self.next_node  = self.Track.GetNode( index + 1 )
        self.this_hit   = self.this_node.hit
        self.prev_state = self.prev_node.filt_state.state
        self.prev_cov   = self.prev_node.filt_state.CovarianceMatrix
        self.MSMatrix   = self.MultipleScatteringMatrix( index )
        self.MMMatrix   = self.MeasurementMatrix( index )
        self.MMMatrixT  = self.MMMatrix.T()
        self.TMatrix    = self.TransportMatrix( index )
        self.TMatrixT   = self.TMatrix.T()
        self.NMatrix    = self.NoiseMatrix( index )
        self.NMatrixI   = self.NMatrix.Inverse()


    def Predict( self, index ):
        '''
            Predict the value for step "index".
        '''
        if not index:
            raise ValueError( 'The initial step cannot be predicted' )

        x_predicted = self.TMatrix ** self.prev_state
        C_predicted = self.TMatrix ** self.prev_cov ** self.TMatrixT + self.MSMatrix
        r_predicted = self.this_hit - self.MMMatrix ** x_predicted
        R_predicted = self.NMatrix + self.MMMatrix ** C_predicted ** self.MMMatrixT
    
        self.this_node.pred_state = KalmanMeasurement( x_predicted, C_predicted )
        self.this_node.pred_resid = KalmanMeasurement( r_predicted, R_predicted )
        
    def Filter( self, index ):
        
        C_filtered = ( self.this_node.pred_state.CovarianceMatrix.Inverse() + self.MMMatrixT ** self.NMatrixI ** self.MMMatrix ).Inverse()
        GainMatrix = C_filtered ** self.MMMatrixT ** NMatrixI
        x_filtered = self.prev_state + GainMatrix ** ( self.this_hit - self.MMMatrix ** self.prev_state )
        
        projector  = Array.Identity( self.ndim ) - self.MMMatrix ** GainMatrix
        r_filtered = projector ** self.this_node.pred_resid
        R_filtered = projector ** self.NMatrix
        
        chi2plus = r_filtered ** R_filtered.Inverse() ** r_filtered
        newchi2  = self.prev_node.cumchi2 + chi2plus
        
        self.this_node.filt_state = KalmanMeasurement( x_filtered, C_filtered )
        self.this_node.filt_resid = KalmanMeasurement( r_filtered, R_filtered )
        self.this_node.chi2       = chi2plus
        self.this_node.cumchi2    = newchi2
    
    def Smooth( self, index ):
        #Still not working
        GainMatrix = self.this_node.filt_state.CovarianceMatrix ** self.TMatrixT ** self.next_node.pred_state.CovarianceMatrix.Inverse()
        x_smooth   = self.this_node.filt_state.state + GainMatrix# ** ( ??? )
    
    def Fit( self ):
        assert self.HaveInitialState and self.HaveMeasurements and self.guess, '''
            You must set the initial state and the measurements. Use the SetInitialState and SetMeasurements methods.
            You must also call the SetInitialGuess method even if you don't have one. If no initial guess is given, it will be estimated by this class.'''
        
        for i in range( 1, self.Track.nnodes ):
            self._NewNode( i )
            self.Predict( i )
            self.Filter( i )
        for i in reversed(range( self.Track.nnodes - 1 )):
            self.NewNode( i )
            self.Smooth( i )

        return self.Track
        





