#!/nfs/hicran/project/compass/analysis/fkrinner/Python_ultra/python_install/bin/python
import sys
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/MassIndependentFit')
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/MassIndependentFit/cards')
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/convertTextOutput')
from create_bootstrap import create_bootstrap
from convertTextOutput import getIntegralMatrixAverage
from PWA import perform_PWA
from shutil import rmtree, copyfile
import numpy as np
import numpy.linalg as la
import os
from write_isobars import update_bootstrap


def get_mixing_indices(cardname):
	with open(cardname,'r') as ininin:
		nwave = 0
		WAS_PREV = False
		old_name_basis = "will_never_come"
		deiso_borders=[]
		max_deiso = 0
		for line in ininin.readlines():
			if line.strip().startswith("*IWAVENAM"):
				nwave+=1
				if "f0_" in line or "f2_" in line or "rho_" in line:
					line_name_basis = line.split("'")[1].split("_")[0]
					if not line_name_basis == old_name_basis:
						deiso_borders.append(nwave)
					old_name_basis = line_name_basis			
					max_deiso = nwave
		deiso_borders.append(max_deiso)
	return deiso_borders

def create_matrix_files(mMin,mMax,binWidth,integrals, target, minindex,maxindex):
	pathfile = target+'/integral_files.dat'
	while len(pathfile)<96:
		pathfile = pathfile+'t' # (*) for some reason, the program does not sccept shorter names
	with open(pathfile,'w') as paths:
		mlow = mMin
		mup  = mMin + binWidth
		while mup <= mMax:
			matrixFile = target+'/matrix_'+str(mlow)+'_'+str(mup)	
			paths.write(str(mlow)+' '+str(mup)+'\n'+matrixFile+'\n')
			matrix_values = get_eigenbasis(mlow,mup,integrals,minindex,maxindex)
			with open(matrixFile,'w') as matrix:
				matrix.write(str(minindex)+' '+str(maxindex)+'\n')	
				for matrix_line in matrix_values:
					for value in matrix_line:
						matrix.write('('+str(value.real)+','+str(value.imag)+') ')
					matrix.write('\n')
			mlow=mup
			mup+=binWidth
	if len(pathfile) < 96: #should be obsolete, due to (*)
		raise IndexError # Length of the matrix file paht must be longer, don't know why

	return pathfile

def create_simple_matrix(mMin,mMax,bin,target,mixing_indices):
	pathfile = target+'/integral_files.dat'
	dim = mixing_indices[-1] - mixing_indices[0] + 1
	matrix_values=[]	
	for i in range(dim):
		matrix_values.append([0.+0.j]*dim)
	for i in range(dim):
		matrix_values[i][i] = 1.+0.j
		if i > 0:
			pass
			matrix_values[i-1][i] = -.5+0.j
			matrix_values[i][i-1] = -.5+0.j
		if i < dim-1:
			pass
			matrix_values[i+1][i] = -.5+0.j
			matrix_values[i][i+1] = -.5+0.j
	for i in range(dim):
#		print i,mixing_indices
		if i + mixing_indices[0] in mixing_indices and not i==0:
			matrix_values[i-1][i] = 0.+0.j
			matrix_values[i][i-1] = 0.+0.j
	matrix_values[dim-1][dim-2] = -.5+0.j
	while len(pathfile)<96:
		pathfile = pathfile+'t' # (*) for some reason, the program does not sccept shorter names	
	with open(pathfile,'w') as paths:
		mlow = mMin
		mup = mMin + bin
		while mup <= mMax:
			matrixFile = target+'/matrix_'+str(mlow)+"_"+str(mup)
			paths.write(str(mlow)+' '+str(mup)+'\n'+matrixFile+'\n')
			with open(matrixFile,'w') as matrix:
				matrix.write(str(mixing_indices[0])+' '+str(mixing_indices[-1])+"\n")
				for matrix_line in matrix_values:
					for value in matrix_line:
						matrix.write('('+str(value.real)+','+str(value.imag)+') ')
					matrix.write("\n")
			mlow=mup
			mup+=bin

	return pathfile
	

def get_eigenbasis(mlow,mup,integrals,minindex,maxindex):
	dim = maxindex-minindex+1 # +1 since the maxindex das to be taken as well
	integral_data = getIntegralMatrixAverage(mlow,mup,integrals,False,False)
	matrix = []
	for i in range(dim):
		matrix.append([0.+0.j]*dim)
	for i in range(minindex,maxindex+1):
		for j in range(minindex,maxindex+1):
			matrix[i-minindex][j-minindex] = integral_data[1][i-1][j-1] # Indices are counted from 1(fortran style), matrix starts at [0][0](Python style)
	small_indices = []
	smalldim = 0
	for i in range(dim):
		if not matrix[i][i] == 0.j:
			small_indices.append(i)
			smalldim+=1
	bigeigenbasis  = []
	for i in range(dim):
		line = [0.+0.j]*dim
		bigeigenbasis.append(line)
	EV_STUFF = False
	if EV_STUFF:
		smallmatrix=[] # Only contains nonzero waves...
		for i in range(smalldim):
			smallmatrixline=[]
			ii = small_indices[i]
			for j in range(smalldim):	
				jj = small_indices[j]
				smallmatrixline.append(matrix[ii][jj])
			smallmatrix.append(smallmatrixline)
		smalleigenbasis = get_small_eigenbasis(smallmatrix)
	
		for i in range(smalldim):
			ii = small_indices[i]
			for j in range(smalldim):
				jj = small_indices[j]
				bigeigenbasis[ii][jj] = smalleigenbasis[i][j]
	else:
		for small_index in small_indices:
			bigeigenbases[small_index][small_index  ] = 1.+0.j
			try:
				
				bigeigenbases[small_index][small_index+1] = -.5+0.j
			except IndexError:
				pass
			try:
				bigeigenbases[small_index][small_index-1] = -.5+0.j
			except IndexError:
				pass
	return bigeigenbasis

def get_small_eigenbasis(matrix,minEV = 0.1):
	dim = len(matrix)
	for i in range(dim):
		for j in range(dim):
			if not i==j:
				matrix[i][j]/=(matrix[i][i]*matrix[j][j])**.5
	for i in range(dim):
		matrix[i][i] = 1.+0.j
	aa = np.zeros((dim,dim),dtype = np.complex128)
	for i in range(dim):
		for j in range(dim):
			aa[i,j] = matrix[i][j]
	val,vec = la.eig(aa)

	vec = np.transpose(vec)
	basis = []
	for i in range(dim):
		basis.append([0.+0.j]*dim)
	for i in range(dim):
		if val[i] > minEV:
			for j in range(dim):
				basis[i][j] = vec[i,j]
		else:
			print "NOTICE:!!!:!!!:!!!:"
			print "One eigenvector omitted:",val[i]
	return basis



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
			ACUTALLY_FIT = True	):
	
	tBins = [tBin]

	card_simple_basis = create_bootstrap( name+'_sim',
				modestring,
				template = template,
				bs_path = path_bootstrap,
				cardfolder = cardfolder)
	
	cardname = card_simple_basis.replace('_sim','')
	mixing_indices = get_mixing_indices(cardfolder+'/'+card_simple_basis)
	matrix_min = mixing_indices[ 0]
	matrix_max = mixing_indices[-1]
	print "Mixing indices are:",matrix_min,matrix_max

	matrix_folder  = target+'/'+name+'/matrix'
	if startStage ==0: # Add creating the matrix files as first stage
#		if ACUTALLY_FIT: # le flage for le debuge
#			perform_PWA(
#				card_simple_basis,
#				name+'_sim',		# Name of the fit
#				mMin,			# Lower mass Limit
#				mMax,			# Upper mass Limit
#				tBins,			# t' bins
#				seeds,			# seeds
#				1,			# 
#				1,			# ONLY INTEGRALS
#				False	,		# 
#				maxResubmit,		# Number of maximum resubmits
#				cleanupWramp,		# Clean
#				cleanupLog,		#	up
#				cleanupFit,		#
#				cleanupInt,		#		flags
#				intBinWidth,		# Bin width for integrals
#				pwaBinWidth,		# Bin width for PWA
#				target,			# Target folder
#				cardfolder,		# Folder with card
#				intSource,		# Source for integrals
#				pwaSource,		# Source for PWA (events)
#				cleanupCore,		# Cleanupt core files
#				MC_Fit,			# Flag if fit to MC events
#				treename,		# Gives the name of the ROOT tree in thr input files Standard name, if none
#				wrampmode,		# Wrampmode, only used, if wramp-files are inclomplete. 
#				COMPENSATE_AMP,		# Compensate_amp flag in the card
#				PRINT_CMD_ONLY		# Flag to print submit commands rather than executing them
#								)
		if not os.path.isdir(target+'/'+name):
			os.makedirs(target+'/'+name)
		if not os.path.isdir(target+'/'+name+'/'+'matrix/'):
			os.makedirs(target+'/'+name+'/'+'matrix/')
#		simIntDir =target+'/'+name+'_sim/integrals/'+tBin[0]+'-'+tBin[1]+'/'
#		pathFile = create_matrix_files(float(mMin),float(mMax),float(pwaBinWidth),simIntDir,matrix_folder,matrix_min, matrix_max)
		pathFile = create_simple_matrix(float(mMin),float(mMax),float(pwaBinWidth),matrix_folder,mixing_indices)
#		if cleanupSim:
#			rmtree(target+'/'+name+"_sim")
		with open(cardfolder+'/'+cardname,'w') as outoutout:
			with open(cardfolder+'/'+card_simple_basis,'r') as ininin:
				for line in ininin.readlines():
					outoutout.write(line)
					if 'PATH_LOAD_BOOTSTRAP' in line:
						outoutout.write('USERDEP_AMP\n')
						outoutout.write("PATH_MIXING_MATRIX '"+pathFile+"'\n")
	if maxStage>0 and proceedStages:
		if ACUTALLY_FIT:
			perform_PWA(
				cardname,
				name,			# Name of the fit
				mMin,			# Lower mass Limit
				mMax,			# Upper mass Limit
				tBins,			# t' bins
				seeds,			# seeds
				max(startStage,1),	# Stage to start at
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
			)
	return target+'/'+name+'/fit/'+tBin[0]+'-'+tBin[1]+'/'


if __name__ == "__main__":
	path_bootstrap = '/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/bootstrap/data_PWA_bs_test'
	bootstrap_result = bootstrap_step(
			modestring		=	'didiiiiii',
			name			=	'PWA_bs_test_6',
			path_bootstrap		=	path_bootstrap,
			template		=	'template_bootstrap_MC.dat',
			mMin			=	'1.50',
			mMax			=	'1.54',
			tBin			=	['0.14077','0.19435'],
			seeds			=	['12345'],
			startStage		=	4,
			maxStage		=	4,
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
			ACUTALLY_FIT		= 	True
	)
	update_bootstrap(bootstrap_result,path_bootstrap, SET_IS_WITH_MATRIX_FALSE=False)


