
from akfbeta import *

alex = Alex('alex')
alex.nevts = 1000

root = ROOTSvc('root',fname='DatosAMirar.root')
alex.addsvc(root) 



mierda = GenerateBeta('myone')
mierda.imports.append('root')

alex.addalg(mierda)
alex.run()




