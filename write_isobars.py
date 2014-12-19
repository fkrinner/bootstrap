import os
import sys
import shutil
sys.path.append('/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/convertTextOutput')
from convertTextOutput import get2D
import numpy as np
import numpy.linalg as la


list_of_possible_key_waves=[
'1-(0-+)0+ rho0mp0pP pi P                                    ',
'1-(1++)0+ rho1pp0pS pi S                                    ',
'1-(1++)1+ rho1pp1pS pi S                                    ',
'1-(2++)1+ rho2pp1pD pi D                                    ',
'1-(2-+)0+ rho2mp0pP pi P                                    ',
'1-(2-+)1+ rho2mp1pP pi P                                    ',
'1-(0-+)0+ f00mp0pS pi S                                     ',
'1-(1++)0+ f01pp0pP pi P                                     ',
'1-(2-+)0+ f02mp0pD pi D                                     ',
'1-(0-+)0+ (pipi)_S pi S                                     ',
'1-(0-+)0+ rho pi P                                          ',
'1-(1++)0+ rho pi S                                          ',
'1-(1++)1+ rho pi S                                          ',
'1-(1++)0+ (pipi)_S pi P                                     ',
'1-(2++)1+ rho pi D                                          ',
'1-(2-+)0+ rho pi P                                          ',
'1-(2-+)1+ rho pi P                                          ',
'1-(2-+)0+ (pipi)_S pi D                                     ']


def getJPC(jpc,data,M='',wave=''):
	"""
	Returns only the part of the data with right JPC
	Both times, the format is the output of 'convertTextOutput.get2D'
	"""
	pnew=[]
	for point in data[0]:
		if point[4] == jpc and (point[13]==M or M=='')and (point[14]==wave or wave==''):
			pnew.append(point[:])
	return [pnew,data[1]]


def getBin(m,binning):
	"""
	Returns the corresponding bin number for m in binning
	"""
	if m < binning[0]:
		raise ValueError("m < m_min")
	for i in range(len(binning)-1):
		if binning[i] <=m and binning[i+1] > m:
			return i
	raise ValueError("m > m_max")

def write_cpp_isobar(data):
	"""
	Returns the data in the format: [binning3[min3],binning2[bin2],re[bin3][bin2],im[bin3][bin2]]
	"""
	binning2 = []
	binning3 = []
	for point in data[0]:
		if not point[0] in binning3:
			binning3.append(point[0])
		if not point[1] in binning3:
			binning3.append(point[1])
		if not point[2] in binning2:
			binning2.append(point[2])
		if not point[3] in binning2:
			binning2.append(point[3])
	binning3.sort()
	binning2.sort()
	line=[0. for i in range(len(binning2)-1)]
	re = [line[:] for i in range(len(binning3)-1)]
	im = [line[:] for i in range(len(binning3)-1)]
	for point in data[0]:	
		m3=(point[0]+point[1])/2
		bin3 = getBin(m3,binning3)
		m2=(point[2]+point[3])/2
		bin2 = getBin(m2,binning2)
		re[bin3][bin2] = point[7]
		im[bin3][bin2] = point[9]
	return[binning3,binning2,re,im]

def write_cpp_lookup_table(filename,data):
	"""
	Writes the data to the files: 
	filename+"_re.dat" for the real part
	filename+"_im.dat" for the imag part
	"""
	try:
		re = data[2]
		im = data[3]
	except:
		print data
		raw_input()
	out_re = open('./'+filename+'_re.dat','w')
	out_im = open('./'+filename+'_im.dat','w')
	out_bin3=open('./'+filename+'_bin3.dat','w')
	out_bin2=open('./'+filename+'_bin2.dat','w')
	for i in range(len(re)):
		for j in range(len(re[i])):
			out_re.write(str(re[i][j])+'   ')
			out_im.write(str(im[i][j])+'   ')
	binning3 = data[0]
	for val in binning3:
		out_bin3.write(str(val)+'   ')
	binning2 = data[1]
	for val in binning2:
		out_bin2.write(str(val)+'   ')
	out_re.close()
	out_im.close()
	out_bin3.close()
	out_bin2.close()

def get_bootstrap_name(jpc,m,isob):
	if jpc == '0-+' and m=='0' and isob == 'rho_':
		return 'rho0mp0pP'
	if jpc == '1++' and m == '0' and isob == 'rho_':
		return 	"rho1pp0pS"	
	if jpc == '2++' and m=='1' and isob == 'rho_':
		return "rho2pp1pD"
	if jpc == '1++' and m =='0' and isob == 'f0_':
		return "f01pp0pP"
	if jpc == '0-+' and m=='0' and isob == 'f0_':
		return "f00mp0pS"
	if jpc == '1++' and m =='1' and isob =='rho_':
		return "rho1pp1pS"
	if jpc == '2-+' and m=='0' and isob == 'rho_':
		return "rho2mp0pP"
	if jpc == '2-+' and m=='1' and isob == 'rho_':
		return "rho2mp1pP"
	if jpc == '2-+' and m=='0' and isob== 'f0_':
		return "f02mp0pD"
	return 'nibp' # Not in bootstrap procedure


def create_bootstrap_files(direct,name):
	"""
	Creates bootsrap files for all de-isobarred JPCs in 'direct'
	"""
	data = get2D(direct)
	jpcs = []
	for point in data[0]:
		if not [point[4],point[13]] in jpcs:
			jpcs.append([point[4],point[13]])
	for jpc in jpcs:
		print "Write tables for JPC = "+jpc[0]+' '+jpc[1]
		write_cpp_lookup_table(name+'_'+jpc[0]+'_'+jpc[1],write_cpp_isobar(getJPC(jpc[0],data,jpc[1],jpc[2])))


def kompleks(strrr):
	reim = strrr.replace('(','').replace(')','').split(',')
	return complex(float(reim[0]),float(reim[1]))


def get_matrix(path_mixing_matrix,m3pi):
		matrixx=[]
		with open(path_mixing_matrix,'r') as matrices: # find the path of the matrix corresponding to the mass bin
			nextis = False
			for line in matrices.readlines():
				if nextis:
					path_matrix = line.strip()
					break
				if float(line.split()[0]) <= m3pi and float(line.split()[1]) > m3pi:
					nextis=True
		with open(path_matrix,'r') as matrix:
			first = True
			for line in matrix.readlines():
				if first:
					first = False
					startindex = int(line.split()[0])
					stoppindex = int(line.split()[1])
				else:
					valline = [kompleks(chunk) for chunk in line.split()]
					matrixx.append(valline)
		return startindex,stoppindex,matrixx

def get_matrix_JPC(sta, sto, ca):
	JPCS = []
	with open(ca,'r') as card:
		nWave = 0
		for line in card.readlines():
			if line.strip().startswith("*IWAVENAM"):
				nWave+=1
				if nWave >= sta:
					namoi = line.split("'")[1]
					jpc = namoi[3:6]
					M   = namoi[7:8]
					if 'f0_' in namoi:
						isob = 'f0_'
					elif 'rho_' in namoi:
						isob = 'rho_'
					else:
						continue	
#						print namoi
#						raise Exception # No isobar 
					JPC = [jpc,M,isob]
					if not JPC in JPCS:
						JPCS.append(JPC)
				if nWave > sto:
					break
	return JPCS

def update_bootstrap(direct,datadir="./data/",SET_IS_WITH_MATRIX_FALSE = False):
	"""
	Updates the bootstrap isobars, taking data from 'direct' and updating the corresponding files in './data/'
	Works only with one mass bin
	"""
	log = open(datadir+'/log','a')
	log.write(direct+'\n')
	print 'get2D('+direct+')'
	for key_wave in list_of_possible_key_waves:
		try:
			print direct
			data = get2D(direct,keywave=key_wave)	
			break
		except:
			print key_wave,"not in the fit"
			pass
	jpcs = []
	IS_WITH_MATRIX = False
	for fn in os.listdir(direct):
		if 'card' in fn:
			le_card = direct+'/'+fn
			with open(le_card,'r') as card:
				for line in card.readlines():
					if "PATH_MIXING_MATRIX" in line:
						path_mixing_matrix = line.split()[1].replace("'","").replace('"','')
						print "Mixing matrix found"
						IS_WITH_MATRIX = True
						break
			break
	if SET_IS_WITH_MATRIX_FALSE:
		IS_WITH_MATRIX = False

	for point in data[0]:
		if not [point[4],point[13],point[14]] in jpcs:
			jpcs.append([point[4],point[13],point[14]])
	m3pi_min = data[0][0][0]
	m3pi_max = data[0][0][1]
	m3pi = (m3pi_min+m3pi_max)/2
	if IS_WITH_MATRIX:
		startindex,stoppindex,matrix = get_matrix(path_mixing_matrix,m3pi)
		JPCS_matrix = get_matrix_JPC(startindex,stoppindex,le_card)
	new_file_names = []
	for jpc in jpcs:
		only_one = write_cpp_isobar(getJPC(jpc[0],data,jpc[1],jpc[2]))
		name = get_bootstrap_name(jpc[0],jpc[1],jpc[2])
		if name == 'nibp':
			print 'Wave not in bootstrap procedure: '+str(jpc)
			continue
		re = only_one[2]
		im = only_one[3]
		out_re = open(datadir+'/'+name+'_re.dat','w')
		out_im = open(datadir+'/'+name+'_im.dat','w')
		new_file_names.append(name)
		for i in range(len(re)):
			intens = 0. #Normalize
			for j in range(len(re[i])):
				intens+=re[i][j]**2.+im[i][j]**2.
			intens**=.5
			for j in range(len(re[i])):
				out_re.write(str(re[i][j]/intens)+'   ')
				out_im.write(str(im[i][j]/intens)+'   ')
		out_re.close()
		out_im.close()

	if IS_WITH_MATRIX:
		complex_values = []
		count_index=0
		index_borders=[0]
		for jpc in JPCS_matrix:
			name = get_bootstrap_name(jpc[0],jpc[1],jpc[2])
			with open(datadir+'/'+name+'_re.dat','r') as reall:
				with open(datadir+'/'+name+'_im.dat','r') as imagg:
					ch_re = reall.read().split()
					ch_im = imagg.read().split()
					if not len(ch_re) == len(ch_im):
						raise IndexError # Number of values does not match
					for i in range(len(ch_re)):
						count_index+=1
						complex_values.append(complex(float(ch_re[i]),float(ch_im[i])))
			index_borders.append(count_index)
		dim = len(complex_values)
		if not len(matrix) == dim:
			raise IndexError # Matrix has wrong dimension
		for i in range(dim):
			if not len(matrix[i]) == dim:
				raise IndexError # matrix not quadratic
		print complex_values
		transformed=[complex(0.,0.)]*len(complex_values)
		for i in range(dim):
			for j in range(dim):
#				matrix = invert(matrix)
				transformed[i] += matrix[j][i]*complex_values[j]
		print "__________________________________________________________________________________________"
		for line in matrix:
			print line

		print "{}{}{}{}{}{}{}{}{}{}{}{}{}"
		print transformed
		for i in range(len(JPCS_matrix)):
			jpc = JPCS_matrix[i]
			name = get_bootstrap_name(jpc[0],jpc[1],jpc[2])
			with open(datadir+'/'+name+'_re.dat','w') as reall:
				with open(datadir+'/'+name+'_im.dat','w') as imagg:
					for j in range(index_borders[i],index_borders[i+1]):
						reall.write(str(transformed[j].real)+' ')
						imagg.write(str(transformed[j].imag)+' ')
	for jpc in jpcs:# copy files to the storage
		name = get_bootstrap_name(jpc[0],jpc[1],jpc[2])
		i=0
		while os.path.isfile(datadir+'/'+name+'_re.dat'+str(i)):
			i+=1
		while os.path.isfile(datadir+'/'+name+'_im.dat'+str(i)):
			i+=1
		shutil.copyfile(datadir+'/'+name+'_re.dat',datadir+'/'+name+'_re.dat'+str(i))
		shutil.copyfile(datadir+'/'+name+'_im.dat',datadir+'/'+name+'_im.dat'+str(i))
		log.write('\t./data/'+name+'_re.dat'+str(i)+'\n')
		log.write('\t./data/'+name+'_im.dat'+str(j)+'\n')
		print name+' bootstrap file written'
	log.close()


def print_bootstrap_waves():
	print "rho1pp0pS"
	print "rho2pp1pD"
	print "f01pp0pP"
	print "f00mp0pS"
	print "rho1pp1pS"
	print "rho0mp0pP"
	print "rho2mp0pP"
	print "rho2mp1pP"
	print "f02mp0pD"









