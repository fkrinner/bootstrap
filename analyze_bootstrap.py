import ROOT
import numpy
import math
from math import pi
from time import sleep
from ROOT import TH1D
from ROOT import TH2D
from ROOT import TGraph
from rootpy.plotting import Canvas
"""
Analyzes the bootstrapping results:
Scipt to be used in interactive mode (ipython)
loads 'n_iter' iterations for all waves in 'waves' (number may be specified separately as well)
- plot_wave(wave): Plots the steps for wave in a TH2D(...)
- list_waves(): Prints all waves in the procedure
- write_wave(wave,mode,suff,one): Writes TH1D(...) for each step of wave. Mode specifies (intensity (default), real, imag, phase, argand), suff specifies the output format, one specifies, whether everything should be in one output file
"""
canv = Canvas(name = "canv0", title = "PWA")

folder = './data_only_1pp0prhoS/'

n_iter = 7

#waves = [['rho0mp0pP',11],['rho1pp0pS',9],['rho1pp1pS',9],['rho2mp0pP',11],['rho2mp1pP',9],['rho2pp1pD',11],['f00mp0pS',5],['f01pp0pP',5],['f02mp0pD',5]]
#waves = [['rho0mp0pP',12],['rho1pp0pS',12],['rho1pp1pS',12],['rho2mp0pP',12],['rho2mp1pP',12],['rho2pp1pD',12],['f00mp0pS',12],['f01pp0pP',12],['f02mp0pD',12]]
waves = [['rho0mp0pP',n_iter],['rho1pp0pS',n_iter],['rho1pp1pS',n_iter],['rho2mp0pP',n_iter],['rho2mp1pP',n_iter],['rho2pp1pD',n_iter],['f00mp0pS',n_iter],['f01pp0pP',n_iter],['f02mp0pD',n_iter]]

f0_waves=['f00mp0pS','f01pp0pP','f02mp0pD']



bin_file = open(folder+'bin_rho.dat','r')
for line in bin_file.readlines():
	chunks=line.split()
	binning_rho = [float(chunk) for chunk in chunks]
bin_file.close()
binning_rho = numpy.asarray(binning_rho,dtype=numpy.float64)		

bin_file = open(folder+'bin_f0.dat','r')
for line in bin_file.readlines():
	chunks=line.split()
	binning_f0 = [float(chunk) for chunk in chunks]
bin_file.close()
binning_f0 = numpy.asarray(binning_f0,dtype=numpy.float64)		

data={}
for wave in waves:
	wave_name = wave[0]
	data[wave_name]=[]
	for i in range(wave[1]):
		dat=[]
		real_file = folder+wave_name+'_re.dat'+str(i)
		imag_file = folder+wave_name+'_im.dat'+str(i)
		real_data = open(real_file,'r')
		print "read: "+real_file
		for line in real_data.readlines():
			chunks=line.split()
			dat.append([float(chunk) for chunk in chunks])
		real_data.close()
		imag_data=open(imag_file,'r')
		print "read: "+imag_file
		for line in imag_data.readlines():
			chunks = line.split()
			dat.append([float(chunk) for chunk in chunks])
		imag_data.close()
		data[wave_name].append(dat)

def list_waves():
	listt=[]
	for wave in waves:
		print wave[0]
		listt.append(wave[0])
	return listt

def plot_wave(wave,mode='int'):
	if wave in f0_waves:
		binning=binning_f0
	else:
		binning=binning_rho
	dat = data[wave]
	nIter = 0
	hist = TH2D('iter_'+str(nIter),'iter_'+str(nIter), len(binning)-1,binning, len(dat),0,len(dat))
	for j in range(len(dat)):
		iteration = dat[j]
		nIter+=1
		for i in range(len(iteration[0])):
			if mode == 'int':
				hist.SetBinContent(i+1,j+1,iteration[0][i]**2+iteration[1][i]**2)
			if mode == 'pha':
				re = iteration[0][i]
				im = iteration[1][i]
				pha = math.atan2(im,re)
				hist.SetBinContent(i+1,j+1,pha)
		if mode == 'pha':
			for i in range(1,len(iteration[0])):
				if hist.GetBinContent(i+1,j+1) == 0.:
					continue
				while hist.GetBinContent(i+1,j+1) > hist.GetBinContent(i,j+1)+pi:
					hist.SetBinContent(i+1,j+1,hist.GetBinContent(i+1,j+1)-2*pi)
				while hist.GetBinContent(i+1,j+1) < hist.GetBinContent(i,j+1)-pi:
					hist.SetBinContent(i+1,j+1,hist.GetBinContent(i+1,j+1)+2*pi)
	hist.Draw()
	return hist

wns=[]
for key in waves:
	wns.append(key[0])
for w in wns:
	print w

def write_wave(wave,mode='int',suff='pdf',one=False):
	if wave in f0_waves:
		binning=binning_f0
	else:
		binning=binning_rho
	dat = data[wave]
	for j in range(len(dat)):
		if mode == 'int' or mode == 'pha':
			hist = TH1D(wave+'_step_'+str(j),wave+'_step_'+str(j),len(binning)-1,binning)
		if mode == 'arg':
			res=[]
			ims=[]
		iteration = dat[j]
		if not mode == 'arg':
			hist.GetXaxis().SetTitle("Mass of the #pi^{#font[122]{-}}#pi^{#font[122]{+}} System (GeV/#it{c}^{2})")
		if mode == 'int':
			hist.GetYaxis().SetTitle('Itensity')
		if mode == 'pha':
			hist.GetYaxis().SetTitle('Phase')
		for i in range(len(iteration[0])):
			if mode == 'int':
				hist.SetBinContent(i+1,iteration[0][i]**2+iteration[1][i]**2)
			if mode == 'pha':
				re = iteration[0][i]
				im = iteration[1][i]
				pha = math.atan2(im,re)
				hist.SetBinContent(i+1,pha)
			if mode == 'arg':
				res.append(iteration[0][i])
				ims.append(iteration[1][i])
		if mode == 'arg':
			res = numpy.asarray(res,dtype=numpy.float64)	
			ims = numpy.asarray(ims,dtype=numpy.float64)
			hist = TGraph(len(res),res,ims)		
		if mode == 'pha':
			for i in range(1,len(iteration[0])):
				if hist.GetBinContent(i+1) == 0.:
					continue
				while hist.GetBinContent(i+1) > hist.GetBinContent(i)+pi:
					hist.SetBinContent(i+1,hist.GetBinContent(i+1)-2*pi)
				while hist.GetBinContent(i+1) < hist.GetBinContent(i)-pi:
					hist.SetBinContent(i+1,hist.GetBinContent(i+1)+2*pi)
		hist.Draw()
		if one:
			if j == len(dat)-1:
				canv.Print('./plots/'+wave+'_'+mode+'.'+suff+')')
			else:	
				canv.Print('./plots/'+wave+'_'+mode+'.'+suff+'(')
		else:
			canv.Print('./plots/'+wave+'_'+mode+'_step_+'+str(j)+'.'+suff)
		canv.Clear()

