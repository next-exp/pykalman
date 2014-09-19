"""
MatPlotLib Implementation of Base class DrawVector
"""

import KDrawVector
import matplotlib.pyplot
import os
import sys
import numpy

#from mpl_toolkits.mplot3d import Axes3D

class MPLDrawVector(KDrawVector.KDrawVector):
    """
    Define an interface to draw hits
    """

    def __init__(self,vector):
        """
        vector =[] or np
        """

        KDrawVector.KDrawVector.__init__(self,vector)
        
        

    def Draw(self):
        """
        Draws Vector
        """
        
        matplotlib.pyplot.show(matplotlib.pyplot.plot(self.vector))
            
        