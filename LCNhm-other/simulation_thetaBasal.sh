MORPHOLOGIES=(sup1 n409 n127 n128 CA196f57-sup cck-28-sup PV-18-sup CA-205R-F57-deep-mEC CA211-LF-59-deep-AMG CA211-RF-57-deep-mPFC cck1-deep PV20-deep PV-63-deep-mEC PV96-deep-mEC);
read -p "Enter Intrinsic (0 1 2 4 6 7 9 12 14 15 16 19 20 21): " INTRINSIC

for ii in {0..1000}; do
	#for INTRINSIC in 0 1 2 4 6 7 9 12 14 15 16 19 20 21; do
	for SYN in 0 1 3; do
		for MORPH in {0..3}; do
			MORPHOLOGY=${MORPHOLOGIES[$MORPH]}
			echo "=========================================="
			echo "       INTRINSIC $INTRINSIC SYN $SYN MORPH $MORPHOLOGY    "
			echo "=========================================="
			sed -i "18s/.*/CELLPROP_MORPHOLOGY = '${MORPHOLOGIES[$MORPH]}'/" LCNhm_configurationfile.py;
			sed -i "19s/.*/CELLPROP_INTRINSIC = $INTRINSIC/" LCNhm_configurationfile.py;
			sed -i "20s/.*/CELLPROP_SYNAPTIC = $SYN/" LCNhm_configurationfile.py;
			python LCNhm_main_argvs.py $INTRINSIC $SYN $MORPHOLOGY;
		done;
	done;
done;

#for ii in {0..1000}; do for INTRINSIC in 0 1 2 4 6 7 9 12 14 15 16 19 20 21; do for SYN in 0 1 3; do echo "=========================================="; echo "       INTRINSIC $INTRINSIC SYN $SYN MORPH $MORPHOLOGY    "; echo "=========================================="; sed -i "18s/.*/CELLPROP_MORPHOLOGY = '$MORPHOLOGY'/" LCNhm_configurationfile.py; sed -i "19s/.*/CELLPROP_INTRINSIC = $INTRINSIC/" LCNhm_configurationfile.py; sed -i "20s/.*/CELLPROP_SYNAPTIC = $SYN/" LCNhm_configurationfile.py; python LCNhm_main_argvs.py $INTRINSIC $SYN $MORPHOLOGY; done; done; done;
