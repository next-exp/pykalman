from ROOT import *
from troot import tcanvas,tconfigure
from troot import tgraph

fname = 'akfnext.root'
t = TFile(fname)
names = ['filter','smooth']
cs = []
tgs = []
for name in names:
	h1,h2 = t.Get('acrev_'+name),t.Get('ac_'+name)
	nn = h1.GetNbinsX()
	x1 = map(lambda i: h1.GetBinContent(i),range(nn))
	x2 = map(lambda i: h2.GetBinContent(i),range(nn))
	a1 = sum(x1)
	a2 = sum(x2)
	xi1 = map(lambda i: sum(x1[:i+1])/a1,range(nn))
	xi2 = map(lambda i: sum(x2[:i+1])/a2,range(nn))
	hs = [h1,h2]
	tg = tgraph([xi1,xi2])
	tconfigure(hs)
	tg.GetXaxis().SetNdivisions(20)
	tg.GetYaxis().SetNdivisions(20)
	cs.append(hs)
	tconfigure(tg)
	cs.append([tg,tg])
	#cs.append(tg)
	#tg.Draw('AL')
	#raw_input('press key to continue ')
c = tcanvas(cs,nx=2,name='ac')
raw_input('press key ')