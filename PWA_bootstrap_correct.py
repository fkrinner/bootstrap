#!/nfs/hicran/project/compass/analysis/fkrinner/Python_ultra/python_install/bin/python
import sys
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/MassIndependentFit')
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/MassIndependentFit/cards')
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/convertTextOutput')
from create_bootstrap import create_bootstrap
from convertTextOutput import getIntegralMatrixAverage
from de_isobarred_corrected import getData
from PWA import perform_PWA
from shutil import rmtree, copyfile
import numpy as np
import numpy.linalg as la
import os
from multiprocessing import Process
import pickle
from random import randint

"""
The name '_correct' does not mean, that its correct, but that is ocrrects the result by the unphsical modes
"""

def bootstrap_step(
			modestring,
			name,
			path_bootstrap,
			template,		# Name of the card
			mMin,			# Lower mass Limit
			mMax,			# Upper mass Limit
			tBin,			# Only one t' bin
			seeds,			# seeds
			startStage,		# Stage to start at
			maxStage,		# Final Stage
			proceedStages,		# Flag to proceed after stages
			maxResubmit,		# Number of maximum resubmits
			cleanupWramp,		# Clean
			cleanupLog,		#	up
			cleanupFit,		#
			cleanupInt,		#		flags
			intBinWidth,		# Bin width for integrals
			pwaBinWidth,		# Bin width for PWA
			target,			# Target folder
			cardfolder,		# Folder with card
			intSource,		# Source for integrals
			pwaSource,		# Source for PWA (events)
			cleanupCore = True,	# Cleanupt core files
			cleanupSim = True,	# Remove the source integrals for the matrix files
			MC_Fit = False,		# Flag if fit to MC events
			treename = None,	# Gives the name of the ROOT tree in thr input files Standard name, if none
			wrampmode = False,	# Wrampmode, only used, if wramp-files are inclomplete. 
			COMPENSATE_AMP	= '0',	# Compensate_amp flag in the card
			PRINT_CMD_ONLY = False,	# Flag to print submit commends rather than executing them
			ACUTALLY_FIT = True,
			replaceMode = 'i'	):
	
	tBins = [tBin]

	cardName1 = create_bootstrap( name,
				modestring,
				template = template,
				bs_path = path_bootstrap,
				cardfolder = cardfolder)

	allCards = [cardName1]
	allNames = [name]

	bindices=[]
	for i in range(len(modestring)):
		if modestring[i] == 'd':
			bindices.append(i)
	
	additional_modestrings=[] 
	if len(bindices) > 1:
		for repbin in bindices:
			addstr = ''
			for i in range(len(modestring)):
				if i == repbin:
					addstr+=modestring[i] # this is a 'd', because this is one bootstrap-mode
				elif i in bindices:
					addstr+=replaceMode
				else:
					addstr+=modestring[i]
			additional_modestrings.append(addstr)
	ccc=0
	for amString in additional_modestrings:
		actName = name+'_'+str(ccc)
		allNames.append(actName)
		allCards.append(create_bootstrap( actName,
				amString,
				template = template,
				bs_path = path_bootstrap,
				cardfolder = cardfolder))
		ccc+=1


	
	if ACUTALLY_FIT:
		procs = []
		
		for i in range(len(allCards)):
			actCard = allCards[i]
			actName = allNames[i]

			proc = Process(target =perform_PWA, args=(
				actCard,
				actName,		# Name of the fit
				mMin,			# Lower mass Limit
				mMax,			# Upper mass Limit
				tBins,			# t' bins
				seeds,			# seeds
				startStage,	# Stage to start at
				maxStage,		# Final Stage
				proceedStages,		# Flag to proceed after stages
				maxResubmit,		# Number of maximum resubmits
				cleanupWramp,		# Clean
				cleanupLog,		#	up
				cleanupFit,		#
				cleanupInt,		#		flags
				intBinWidth,		# Bin width for integrals
				pwaBinWidth,		# Bin width for PWA
				target,			# Target folder
				cardfolder,		# Folder with card
				intSource,		# Source for integrals
				pwaSource,		# Source for PWA (events)
				cleanupCore,		# Cleanupt core files
				MC_Fit,			# Flag if fit to MC events
				treename,		# Gives the name of the ROOT tree in thr input files Standard name, if none
				wrampmode,		# Wrampmode, only used, if wramp-files are inclomplete. 
				COMPENSATE_AMP,		# Compensate_amp flag in the card
				PRINT_CMD_ONLY,		# Flag to print submit commends rather than executing them
			))
			proc.start()
			procs.append(proc)
		for proc in procs:
			proc.join()

	results = []
	for name in allNames:
		results.append(target+'/'+name+'/fit/'+tBin[0]+'-'+tBin[1]+'/')
	return results

def find_text_fits(path):
	"""Gets the folder with textoutput results. Assumes only one folder"""
	folder = None
	for fn in os.listdir(path):
		if os.path.isdir(path+'/'+fn) and fn.startswith('text_fit_'):
			if not folder:
				folder = fn
			else:
				raise IOError # More than one folder found
	return folder

def get_bootstrap_name(key):
	"""converts from usual notation to bootstrap notation"""
	jpc = key[3:6]
	m = key[7]
	isob = key.split()[-1]
	if jpc == '0-+' and m=='0' and isob == 'rho':
		return 'rho0mp0pP'
	if jpc == '1++' and m == '0' and isob == 'rho':
		return 	"rho1pp0pS"	
	if jpc == '2++' and m=='1' and isob == 'rho':
		return "rho2pp1pD"
	if jpc == '1++' and m =='0' and isob == 'f0':
		return "f01pp0pP"
	if jpc == '0-+' and m=='0' and isob == 'f0':
		return "f00mp0pS"
	if jpc == '1++' and m =='1' and isob =='rho':
		return "rho1pp1pS"
	if jpc == '2-+' and m=='0' and isob == 'rho':
		return "rho2mp0pP"
	if jpc == '2-+' and m=='1' and isob == 'rho':
		return "rho2mp1pP"
	if jpc == '2-+' and m=='0' and isob== 'f0':
		return "f02mp0pD"
	return 'nibp' # Not in bootstrap procedure


def updateBootstrap(path_bootstrap,bs_name,data):
	file_name_re = path_bootstrap+'/'+bs_name+'_re.dat'
	file_name_im = path_bootstrap+'/'+bs_name+'_im.dat'
	number_re = 0
	while True:
		store_fn_re = file_name_re.replace('_re.dat','_re'+str(number_re)+'.dat')
		if not os.path.isfile(store_fn_re):
			copyfile(file_name_re,store_fn_re)
			break
		else:
			number_re+=1
	number_im = 0
	while True:
		store_fn_im = file_name_im.replace('_im.dat','_im'+str(number_im)+'.dat')
		if not os.path.isfile(store_fn_im):
			copyfile(file_name_im,store_fn_im)
			break
		else:
			number_im+=1
	with open(path_bootstrap+'/log.log','a') as log:
		log.write("Updated: "+bs_name+", old amplitudes stored in:\n"+store_fn_re+'\n'+store_fn_im+'\n')
	with open(file_name_re,'w') as out_re:
		with open(file_name_im,'w') as out_im:
			for val in data:
				out_re.write(str(val.real)+' ')
				out_im.write(str(val.imag)+' ')





def step_bootstrap(name,path_bootstrap,modestring):

	bootstrap_result = bootstrap_step(
			modestring		=	modestring,
			name			=	name,
			path_bootstrap		=	path_bootstrap,
			template		=	'template_bootstrap_MC.dat',
			mMin			=	'1.50',
			mMax			=	'1.54',
			tBin			=	['0.14077','0.19435'],
			seeds			=	[str(randint(1,100000)) for i in range(10)],
			startStage		=	1,
			maxStage		=	5,
			proceedStages		=	True,
			maxResubmit		=	10,
			cleanupWramp		=	True,
			cleanupLog		=	True,
			cleanupFit		=	True,
			cleanupInt		=	True,
			intBinWidth		=	'0.010',
			pwaBinWidth		=	'0.040',
			target			=	'/nfs/mds/user/fkrinner/massIndepententFits/fits/',
			cardfolder		=	'/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/MassIndependentFit/cards',
			intSource		=	'/nfs/nas/data/compass/hadron/2008/comSkim/MC/PS-MC/trees_for_integrals/m-bins/0.100-1.000/',
			pwaSource		=	'/nfs/nas/data/compass/hadron/2008/comSkim/2008-binned/all/skim_2012',
			cleanupCore 		= 	True,
			cleanupSim		=	False,
			MC_Fit 			= 	False,
			treename 		= 	None,
			wrampmode 		= 	False,	
			COMPENSATE_AMP		= 	'0',
			PRINT_CMD_ONLY		= 	False,
			ACUTALLY_FIT		= 	True,
			replaceMode		=	'b'
	)
	for i in range(len(bootstrap_result)):
		bootstrap_result[i]+='/'+find_text_fits(bootstrap_result[i])

	correctedData = getData(bootstrap_result[0],bootstrap_result[1:],WRITE = False)
#	with open('pickled','w') as outoutout:
#		pickle.dump(correctedData,outoutout)
#	print
	for key in correctedData.iterkeys():
		updateBootstrap(path_bootstrap,get_bootstrap_name(key),correctedData[key])

#	update_bootstrap(bootstrap_result,path_bootstrap, SET_IS_WITH_MATRIX_FALSE=False)

if __name__ == "__main__":
	release_order = [
				"bbddbbbbd",
				"ddbbdbbbb",
				"bbbbbdddb",
				"dbdbbdbbb",
				"bdbdbbdbb",
				"bbbbdbbdd",
				"bddbbbbdb",
				"bbbdddbbb",
				"dbbbbbdbd",
				"bbdbdbdbb",
				"dbbdbbbdb",
				"bdbbbdbbd",
				"bbddbbbbd",
				"ddbbdbbbb",
				"bbbbbdddb",
				"dbdbbdbbb",
				"bdbdbbdbb",
				"bbbbdbbdd",
				"bddbbbbdb",
				"bbbdddbbb",
				"dbbbbbdbd",
				"bbdbdbdbb",
				"dbbdbbbdb",
				"bdbbbdbbd"]

#	name = 'PWA_bs_test_6'
	name = 'PWA_bs_test_7'
	path_bootstrap = '/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/bootstrap/data_PWA_bs_test_7'
	i=0
	for step in release_order:
		step_bootstrap(name+'_'+str(i),path_bootstrap,step)
		i+=1
