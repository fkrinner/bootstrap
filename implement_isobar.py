import os

"""
Implements a new isobar in the COMPASSPWA fortran code
"""

cbw_isob="/nfs/hicran/project/compass/analysis/fkrinner/compassPWA/phys/3pic/compass/2008florian3/cbw_isob.f"
calc_ampl="/nfs/hicran/project/compass/analysis/fkrinner/compassPWA/phys/3pic/compass/2008florian3/calc_ampl_3body_define_3pic.f"

def update_files():
	if os.path.isfile(cbw_isob) and os.path.isfile(cbw_isob+'_new'):
		os.rename(cbw_isob,cbw_isob+'_old')
		os.rename(cbw_isob+'_new',cbw_isob)
	if os.path.isfile(calc_ampl) and os.path.isfile(calc_ampl+'_new'):
		os.rename(calc_ampl,calc_ampl+'_old')
		os.rename(calc_ampl+'_new',calc_ampl)


def is_commented(line):
	if line.startswith('c') or line.startswith('C'):
		return True
	return False
		


def write_cbw(number, name, spin, inter = 5.):
	cbw_in = open(cbw_isob,'r')
	cbw_out= open(cbw_isob+'_new','w')
	for line in cbw_in.readlines():
		if 'MARKER_FOR_PYTHON_SCRIPT' in line:
			cbw_out.write('\n')
			cbw_out.write('      nam_isob('+str(number)+")  = '"+name+"'\n")
			cbw_out.write('      ks_isob('+str(number)+')  = '+str(spin)+'\n')
			cbw_out.write('      rint_isob('+str(number)+') = '+str(inter)+'\n')
			cbw_out.write('\n')
		cbw_out.write(line)
	cbw_in.close()
	cbw_out.close()

def remove_cbw(number):
	cbw_in = open(cbw_isob,'r')
	cbw_out= open(cbw_isob+'_new','w')
	for line in cbw_in.readlines():
		erase_line = False
		if 'nam_isob' in line or 'ks_isob' in line or 'rint_isob' in line:
			chunks = line.split('(')
			if len(chunks)>1:
				if chunks[1].strip().startswith(str(number)):
					erase_line=True
		if not erase_line:
			cbw_out.write(line)
	cbw_in.close()
	cbw_out.close()


def write_ampl(number, spin, name=''):
	in_rho = False
	in_f   = False
	calc_in = open(calc_ampl,'r')
	calc_out= open(calc_ampl+'_new','w')
	to_write=[]
	for line in calc_in.readlines():
		if 'PYTHON_STOP_RHO' in line:
			in_rho=False
			n_bef=0
			for to_write_line in to_write:
				if 'icode_reson' in to_write_line and not is_commented(to_write_line):
					n_bef+=1
			if spin ==1:
				to_write.append('        icode_reson('+str(n_bef+1)+',1,1) = '+str(number))
				if not name == '':
					to_write.append('     ! '+name)
				to_write.append('\n')
			for to_write_line in to_write:
				if 'nreson' in to_write_line and spin ==1:
					calc_out.write('        nreson(  1,1) = '+str(n_bef+1)+'\n')
				else:
					calc_out.write(to_write_line)
			to_write=[]
		if 'PYTHON_STOP_F' in line:
			in_f = False
			n_bef=0
			for to_write_line in to_write:
				if 'icode_reson' in to_write_line and not is_commented(to_write_line):
					n_bef+=1
			if spin == 0 or spin ==2:
				to_write.append('        icode_reson('+str(n_bef+1)+',1,2) = '+str(number))
				if not name == '':
					to_write.append('     ! '+name)
				to_write.append('\n')
			for to_write_line in to_write:
				if 'nreson' in to_write_line and (spin == 0 or spin ==2):
					calc_out.write('        nreson(  1,2) = '+str(n_bef+1)+'\n')
				else:
					calc_out.write(to_write_line)
			to_write=[]
		if not in_f and not in_rho:
			calc_out.write(line)
		else:
			to_write.append(line)
		if 'PYTHON_START_RHO' in line:
			in_rho=True
		if 'PYTHON_START_F' in line:
			in_f=True
	calc_in.close()
	calc_out.close()

def remove_ampl(number):
	in_python = False
	calc_in = open(calc_ampl,'r')
	sect = '0'
	for line in calc_in.readlines():
		if 'PYTHON_START' in line:
			in_python = True
		if 'PYTHON_STOP' in line:
			in_python = False
		if in_python:
			if 'icode_reson' in line and not is_commented(line):
				chunks = line.split('=')
				if len(chunks) > 1:
					if chunks[1].strip().startswith(str(number)):
						sect = chunks[0].strip()[-2]
	if sect == '0':
		print "Section could not be found"
		return
	calc_in.close()
	calc_in = open(calc_ampl,'r')
	calc_out= open(calc_ampl+'_new','w')
	in_python = False
	isob_count=0
	for line in calc_in.readlines():
		write = True
		if 'PYTHON_START' in line:
			in_python = True
		if 'PYTHON_STOP' in line:
			in_python = False
		if in_python:
			if 'nreson' in line and not is_commented(line):
				sect_act = line.split('=')[0].strip()[-2]
				if sect_act == sect:
					nact = int(line.split('=')[1].strip().split()[0])
					calc_out.write('        nreson(  1,'+sect+') = '+str(nact-1)+'\n')
					write= False
			if 'icode_reson' in line and not is_commented(line):
				chunks = line.split('=')
				if len(chunks)>1:
					if chunks[1].strip().startswith(str(number)):	
						write=False
				sect_act = line.split('=')[0].strip()[-2]
				if sect_act==sect and write:
					isob_count+=1
					calc_out.write('        icode_reson('+str(isob_count)+',1,'+line.split(',')[2])
					write=False
		if write:
			calc_out.write(line)

	calc_in.close()
	calc_out.close()


				

def implement_isobar(number, name, spin, inter=5.):
	if number < 3000 or number > 4000:
		print "Isobar number not in the allowed range."
		return
	write_cbw(number, name, spin, inter)
	write_ampl(number, spin, name)
	update_files()


def remove_isobar(number):
	remove_cbw(number)
	remove_ampl(number)
	update_files()




