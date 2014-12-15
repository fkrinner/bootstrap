# Bootstrap script, that does three steps at once, uses fixed names and does some checks

	# For future reference:
	#  the 'steps' just contain names for the sigle stept, used to identify the fit
	#  in principle these names could be anything
	#
	#  the bs_strings contain infomation about the actions to be performed in the according step:
	#  the single chars apply to the wave on the according position. Their meanings are:
	#  i: Usual isobar
	#  b: bootstrapped isobar
	#  d: de-isobarred
	#  0: exclude wave
	#  The different positions belong to the following waves:
	#  rho1pp0pS 	
	#  rho2pp1pD
	#  f01pp0pP
	#  f00mp0pS
	#  rho1pp1pS
	#  rho0mp0pP
	#  rho2mp0pP
	#  rho2mp1pP
	#  f02mp0pD
	#  For these steps, the following actions are performed:
	#  create_bootstrap.py $step $bs_string creates the necessary card and addwave file
	#  
	#  call_perform_pwa.py $step performs the actual PWA
	#  
	#  update_bootstrap.py /nfs/mds/user/fkrinner/massIndepententFits/fits/$step/fit/0.14077-0.19435 loads the de-isobarred results from the fit into the bootstrap files
	#  
nErr=0




if [ $# -gt 1 ]; then
	name_base=$1
	card_template=$2
	if [ $# -gt 2 ]; then # defines a list of release orders
		release_order=$3
	else
		release_order=./standard_release_order.sh
	fi
	source $release_order
	if [ $# -gt 3 ]; then
		n_step=$4 # 3*nStep fits will be performed, since three are always released simultaniously
	else
		n_string=${#bs_strings[@]}
		((n_step=n_string/3))
		echo $n_step steps will be performed
	fi
	if [ $# -gt 4 ]; then
		start_step=$5
	else
		start_step=0
	fi
	datadir=/nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/bootstrap/data_$name_base
	if [[ ! -d $datadir ]]; then
		mkdir -v $datadir 
		source ./copy_start_isobars.sh $datadir
	else
		echo data folder already exists
		echo Bootstrapping script has ended unsuccessfully.
		return
	fi

	cp -v /nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/MassIndependentFit/cards/$card_template $datadir
	for (( l = $start_step ; l < $n_step  ; l++ )); do
		((l1 =  3 * l ))
		((l2 = l1 + 1 ))
		((l3 = l1 + 2 ))
		step1=$name_base"_"$l1
		step2=$name_base"_"$l2
		step3=$name_base"_"$l3
		bs_string1=${bs_strings[$l1]}
		bs_string2=${bs_strings[$l2]}
		bs_string3=${bs_strings[$l3]}
		if [ ${#bs_strings[@]} -le $l3 ]; then
			echo Not enough bs_strings provided, abort.
			((nErr+=1))
			break
		fi
		if [[ -d /nfs/mds/user/fkrinner/massIndepententFits/fits/$step1/ || -d /nfs/mds/user/fkrinner/massIndepententFits/fits/$step2/ || -d /nfs/mds/user/fkrinner/massIndepententFits/fits/$step3/ ]]; then
			echo Target directory exists. Abort procedure.
			((nErr+=1))
			break
		fi
####################################################################################################################################
		source ./mark_pwa.sh $step1 $bs_string1 $card_template $datadir &
		source ./mark_pwa.sh $step2 $bs_string2 $card_template $datadir &
		source ./mark_pwa.sh $step3 $bs_string3 $card_template $datadir &
		while ( [ ! -f ./markers/$step1  ] || [ ! -f ./markers/$step2  ] || [ ! -f ./markers/$step3  ] )
		do
			echo waiting...
			sleep 300
		done
		rm -f ./markers/$step1
		rm -f ./markers/$step2
		rm -f ./markers/$step3
		python update_bootstrap.py /nfs/mds/user/fkrinner/massIndepententFits/fits/$step1/fit/0.14077-0.19435 $datadir # update the bootrap isobars
		python update_bootstrap.py /nfs/mds/user/fkrinner/massIndepententFits/fits/$step2/fit/0.14077-0.19435 $datadir
		python update_bootstrap.py /nfs/mds/user/fkrinner/massIndepententFits/fits/$step3/fit/0.14077-0.19435 $datadir
####################################################################################################################################
	done
	step_fin=$name_base"_fin"
	source ./mark_pwa.sh $step_fin bbbbbbbbb $card_template $datadir  # Change here, if number of bs-isobars changes
	rm -f ./markers/$step_fin
else
	echo No arguments given.
	((nErr+=1))
fi

rm -f bestFits.txt
cd /nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/bootstrap



rm -r /nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/MassIndependentFit/cards/*$name_base.dat # remove all addwaves and cards


if [ $nErr -eq 0 ]; then
	echo Bootstrapping script has ended successfully.
else
	echo Bootstrapping script has ended unsuccessfully.
fi

