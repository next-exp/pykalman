class KFTrack:
    '''
        A simple class to contain all the relevant information relative to the Kalman Filter data.
    '''
    def __init__( self, isbb = False ):
        self.isbb     = bool(isbb)
        self.vertex   = None
        self.zstatesf = []
        self.zstatesb = []
        self.tstatesf = []
        self.tstatesb = []
        self.fstatesf = []
        self.fstatesb = []
        self.sstatesf = []
        self.sstatesb = []
        self.fchif    = []
        self.fchib    = []
        self.schif    = []
        self.schib    = []


from alex import IAlg
import shelve

class KFTrackFiller(IAlg):
    '''
        Filling KFTrack algorithm.
    '''
    def __init__( self, finput, binput, output, isbb ):
        IAlg.__init__( self, input )
        self.finput = finput
        self.binput = binput
        self.output = output
        self.isbb   = isbb
    
    def initialize( self ):
        self.outputfile = shelve.open( self.output )
        self.counter = 0
    
    def execute( self ):
        kftrack = KFTrack( self.isbb )
        
        kfdata = self.evt[self.finput]
        if (not kfdata ): return False
        else: kfdata = kfdata[0]
        for node in kfdata.nodes:
            true = node.getstate('true')
            filt = node.getstate('filter')
            smoo = node.getstate('smooth')
            
            kftrack.zstatesf.append( true.zrun )
            kftrack.tstatesf.append( [i for i in true.vec] )
            kftrack.fstatesf.append( [i for i in filt.vec] )
            kftrack.sstatesf.append( [i for i in smoo.vec] )
            kftrack.fchif   .append( node.getchi2('filter') )
            kftrack.schif   .append( node.getchi2('smooth') )

        kfdata = self.evt[self.binput]
        if (not kfdata ): return False
        else: kfdata = kfdata[0]
        for node in kfdata.nodes:
            true = node.getstate('true')
            filt = node.getstate('filter')
            smoo = node.getstate('smooth')
            
            kftrack.zstatesb.append( true.zrun )
            kftrack.tstatesb.append( [i for i in true.vec] )
            kftrack.fstatesb.append( [i for i in filt.vec] )
            kftrack.sstatesb.append( [i for i in smoo.vec] )
            kftrack.fchib   .append( node.getchi2('filter') )
            kftrack.schib   .append( node.getchi2('smooth') )

        self.outputfile[ str(self.counter) ] = kftrack
        self.counter += 1
        
    def finalize( self ):
        self.outputfile['N'] = self.counter
        self.outputfile.close()
