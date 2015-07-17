# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 19:13:51 2014

@author: brais
"""

import numpy as np
import ROOT as rt
from eficiencia import *


names = ['Chi2MeanFor',
         'Chi2MeanRev',
         'norIn',
         'revIn',
         'asimI',
         'norMed',
         'revMed',
         'asimM',
         'norFin',
         'revFin',
         'asimF',
         'asimasim',
         'asimNorm',
         'asimRev',
         'asimBIn',
         'asimBFin',
         'asimBMed']
lims = [[0,15],
        [0,15],
        [0,15],
        [0,15],
        [-1.05,1.05],
        [0,15],
        [0,15],
        [-1.05,1.05],
        [0,15],
        [0,15],
        [-1.05,1.05],
        [-1.05,1.05],
        [-1.05,1.05],
        [-1.05,1.05],
        [-1.05,1.05],
        [-1.05,1.05],
        [-1.05,1.05]]
        
histf = []
histp = []
canvas = []

for i in range(len(names)):
   histf.append(rt.TH1D(names[i]+'F',names[i]+'F',100,lims[i][0],lims[i][1]))  
   histp.append(rt.TH1D(names[i]+'P',names[i]+'P',100,lims[i][0],lims[i][1])) 
         

infile1 = 'Single_Irene_15_0T_10.dat'                            # Leo los datos del fichero
infile2 = 'Double_Irene_15_0T_10.dat' 



datf = np.loadtxt(infile1,dtype=float)[:]
datp = np.loadtxt(infile2,dtype=float)[:]

Nf,mf = datf.shape
Np,mp = datp.shape

rand = rt.TRandom3(0)

datosMedf = []
datosMedp = []
datosEstf = []
datosEstp = []

for i in range(Np):
    if not datp[i,0]>7: 
        for j in range(mp):
            histp[j].Fill(datp[i,j])
        if rand.Rndm()<=0.5:
            datosMedp.append(datp[i,:])
        else:
            datosEstp.append(datp[i,:])
       
for i in range(Nf):
    if not datf[i,0]>7: 
        for j in range(mf):
            histf[j].Fill(datf[i,j])
        if rand.Rndm()<=0.5:
            datosMedf.append(datf[i,:])
        else:
            datosEstf.append(datf[i,:])

datosMedp = np.array(datosMedp)
datosMedf = np.array(datosMedf)

datosEstp = np.array(datosEstp)
datosEstf = np.array(datosEstf)


Nfm = datosMedf.shape[0]
Npm = datosMedp.shape[0]

meanf = []
meanp = []
covf = np.zeros([mf,mf])
covp = np.zeros([mf,mf])


for i in range(mf):
    xd = datosMedf[:,i]
    xmean = np.mean(xd)
    meanf.append(xmean)    
    
    for j in range(mf):
        yd = datosMedf[:,j]
        ymean = np.mean(yd)
        covf[i,j] = np.mean((yd-ymean)*(xd-xmean))
        
for i in range(mp):
    xd = datosMedp[:,i]
    xmean = np.mean(xd)
    meanp.append(xmean)    
    
    for j in range(mf):
        yd = datosMedp[:,j]
        ymean = np.mean(yd)
        covp[i,j] = np.mean((yd-ymean)*(xd-xmean))
        
meanf = np.array(meanf)
meanp = np.array(meanp)

U = np.dot(np.linalg.inv(covf+covp),(meanf-meanp))
hist1 = rt.TH1D('F','F',100,-20,20)
hist2 = rt.TH1D('F','F',100,-20,20)

cff = 0
cfp = 0
cpf = 0
cpp = 0

Fref = -10.7
F1 = []
F2 = []
for i in range(len(datosEstf)):
    F =np.dot(U,datosEstf[i])
    F1.append(F)
    if F > Fref:
        cff += 1
    else:
        cfp +=1
    hist1.Fill(F)

for i in range(len(datosEstp)):
    F =np.dot(U,datosEstp[i])
    F2.append(F)
    if F > Fref:
        cpf += 1
    else:
        cpp +=1
    hist2.Fill(F)  
cF = rt.TCanvas()
hist1.Draw()
hist2.SetLineColor(2)
hist2.Draw('SAME')



inst = EffPlots(F2,F1)
a = inst.Do(0.05)

outfile = rt.TFile("Fisher10sep_Irene_15_0T.root","recreate")

for i in range(len(names)):
    c = rt.TCanvas()
    
    histf[i].Draw()
    histp[i].Draw('SAME')
    
    histp[i].SetLineColor(2)
    
    canvas.append(c)
    c.Write()
    
cF.Write()
a.Write()
outfile.Close()

'''
print 'Tenemos ',cff,' eventos de la poblacion f identificados correctamente','\n'
print 'Tenemos ',cfp,' eventos de la poblacion f identificados como de la poblacion p','\n'
print 'Tenemos ',cpf,' eventos de la poblacion p identificados como de la poblacion f','\n'
print 'Tenemos ',cpp,' eventos de la poblacion p identificados correctamente','\n'

print 'Por tanto tenemos un porcentaje de exito del ',round(cff*100./Nf,2),'% en la poblacion f y de un ', round(cpp*100./Np,2),'% para la poblacion p'
'''
