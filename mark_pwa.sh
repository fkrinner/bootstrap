#Performs bootstrap PWA woth name $1, bootstrap specification $2 and template card $3. When the PWA is finished, it creates a marker in './markers/' with name $1, so it can be checked
step=$1
bs_string=$2
card_template=$3
datadir=$4
cd /nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/MassIndependentFit/cards # go to cardfolder
python create_bootstrap.py $step $bs_string $card_template $datadir # create the card according to to bs_string (the i's belong to f0 waves, which are left out here, d: de-isobarred b: bootstrapped
cd /nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/MassIndependentFit # got to PWA folder
python call_perform_pwa.py $step # perfom PWA
cd /nfs/hicran/project/compass/analysis/fkrinner/fkrinner/trunk/massDependentFit/scripts/bootstrap # got back to bootstrap folder
sleep `echo $RANDOM%20 | bc`
touch ./markers/$step
