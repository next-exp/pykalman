'''
    This code read Josh's tracks and performs the fit with the KalmanFilter saving all the relevant information. It does not fit the alex philosophy because of the Josh's structure.
'''


from trackdefs import *
from KFBase import KFVector
from ToyParticle import ToyParticle
from KFWolinFilter import KFWolinFilter
from KTrackFitter import KTrackFitter

from alex import Alex
from ialex import IAlg

class JJKalmanFilter( IAlg ):
    
    def __init__( self, measurements = [], name = 'JJKalmanFilter' ):
        IAlg.__init__( self, name )

    def define( self ):
        self.fnb_trk = "{0}/{1}".format(trk_outdir,trk_name)
        self.fnb_fit = "{0}/{1}".format(fit_outdir,trk_name)
        self.files = [ "{0}/{1}_{2}.dat".format(self.fnb_trk,trk_name,ntrk) for ntrk in range(1) ]
    
    def initialize( self ):
        pass
    
    def execute( self ):
        for ntrk,file in enumerate(self.files):
            tpart = ToyParticle( file, False, np.array([sigma_xm,sigma_ym,0.0]), 0 )
            # Set up a KTrackFitter to fit this track, and fit the requested number of times.
#            fSegments = []
            fHits = tpart.SmearedHits(False)
            for ft in range(nfits):
                
#                print "--> Fit {0} of {1}".format(ft,nfits)
                
                # Create a KalmanFilter to be used to do the fitting
 #               logging.debug("-- Creating KFWolinFilter...")
                kfilter = KFWolinFilter("KFWolinFilter");
                
#                logging.debug("-- Creating KTrackFitter...")
                tfitter = KTrackFitter(kfilter,fHits,np.array([sigma_xm,sigma_ym]),tpart,chi2_lim,10.)
                
                # Perform the fit.
#                logging.debug("-- Performing fit...")
                fHits = tfitter.Fit()
                ftrack = tfitter.Track
                print ftrack
                self.evt[self.name] = ftrack
                
                return True
                # Save the segments.
#                fSegments = tfitter.Segments
#            
#            # Write the fit files.
#            f_ftrk = open("{0}/fit_{1}_{2}.dat".format(self.fnb_fit,trk_name,ntrk),"w")
#            f_fseg = open("{0}/seg_{1}_{2}.dat".format(self.fnb_fit,trk_name,ntrk),"w")
#            f_ftrk.write("# segID k x0 y0 z0 p1p p2p p3p p4p chi2p p1f p2f p3f p4f chi2f cfxy cftxy\n")
#            f_fseg.write("# segID nPts chi2avg chi2min chi2max\n")
#            
#            segments = fSegments
#            
#            for seg in segments:
#                
#                # Get the mean, min, and max chi2 values for each segment.
#                #  Set to -1 if the segment is not >= 3 points.
#                mean_chi2 = -1.; min_chi2 = -1.; max_chi2 = -1.;
#                if(len(seg.seg_k) >= 2):
#                    mean_chi2 = np.mean(seg.seg_fchisq[1:]);
#                    min_chi2 = min(seg.seg_fchisq[1:]);
#                    max_chi2 = max(seg.seg_fchisq[1:]);
#                
#                # Print the segment file line.
#                f_fseg.write("{0} {1} {2} {3} {4}\n".format(seg.seg_id,len(seg.seg_k),mean_chi2,min_chi2,max_chi2));
#                
#                # Print the track file lines for each segment point: include the smeared hits.
#                for k,x0,y0,z0,p1p,p2p,p3p,p4p,chi2p,p1f,p2f,p3f,p4f,chi2f,cfxy,cftxy,edep in zip(seg.seg_k,seg.seg_x0,seg.seg_y0,seg.seg_z0,seg.seg_p1p,seg.seg_p2p,seg.seg_p3p,seg.seg_p4p,seg.seg_pchisq,seg.seg_p1f,seg.seg_p2f,seg.seg_p3f,seg.seg_p4f,seg.seg_fchisq,seg.seg_cfxy,seg.seg_cftxy,seg.seg_edep):
#                    if(isinstance(chi2f,KFVector)):
#                        chi2 = chi2f[2];
#                    else:
#                        chi2 = chi2f;
#                    f_ftrk.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17}\n".format(seg.seg_id,k,x0,y0,z0,p1p,p2p,p3p,p4p,chi2p,p1f,p2f,p3f,p4f,chi2f,cfxy,cftxy,edep));
#        return True

    def finalize( self ):
        return

alex = Alex('alex')
kf   = JJKalmanFilter()
alex.addalg(kf)
alex.nevts = 1
alex.run()
