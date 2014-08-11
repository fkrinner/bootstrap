if [[ $# -gt 0 ]]; then
	DATAFOLDER=$1
else
	DATAFOLDER=./data
fi


cp -v ./start_isobars/bin* $DATAFOLDER
echo -------------------------------
cp -v ./start_isobars/wrong_f0_re.dat  $DATAFOLDER/f00mp0pS_re.dat
cp -v ./start_isobars/wrong_f0_im.dat  $DATAFOLDER/f00mp0pS_im.dat
echo -------------------------------
cp -v ./start_isobars/wrong_sig_re.dat $DATAFOLDER/f01pp0pP_re.dat
cp -v ./start_isobars/wrong_sig_im.dat $DATAFOLDER/f01pp0pP_im.dat
echo -------------------------------
cp -v ./start_isobars/wrong_sig_re.dat $DATAFOLDER/f02mp0pD_re.dat
cp -v ./start_isobars/wrong_sig_im.dat $DATAFOLDER/f02mp0pD_im.dat
echo -------------------------------
cp -v ./start_isobars/wrong_rho_re.dat $DATAFOLDER/rho0mp0pP_re.dat
cp -v ./start_isobars/wrong_rho_im.dat $DATAFOLDER/rho0mp0pP_im.dat
echo -------------------------------
cp -v ./start_isobars/wrong_rho_re.dat $DATAFOLDER/rho1pp0pS_re.dat
cp -v ./start_isobars/wrong_rho_im.dat $DATAFOLDER/rho1pp0pS_im.dat
echo -------------------------------
cp -v ./start_isobars/wrong_rho_re.dat $DATAFOLDER/rho1pp1pS_re.dat
cp -v ./start_isobars/wrong_rho_im.dat $DATAFOLDER/rho1pp1pS_im.dat
echo -------------------------------
cp -v ./start_isobars/wrong_rho_re.dat $DATAFOLDER/rho2mp0pP_re.dat
cp -v ./start_isobars/wrong_rho_im.dat $DATAFOLDER/rho2mp0pP_im.dat
echo -------------------------------
cp -v ./start_isobars/wrong_rho_re.dat $DATAFOLDER/rho2mp1pP_re.dat
cp -v ./start_isobars/wrong_rho_im.dat $DATAFOLDER/rho2mp1pP_im.dat
echo -------------------------------
cp -v ./start_isobars/wrong_rho_re.dat $DATAFOLDER/rho2pp1pD_re.dat
cp -v ./start_isobars/wrong_rho_im.dat $DATAFOLDER/rho2pp1pD_im.dat
echo -------------------------------
