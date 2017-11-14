#!/bin/bash

mkdir "$TMPDIR"/VARsim_UT/

cd "$HOME"/VARsim_UT

for i in `seq 1 100`;
do
	sed s/iter/$i/g VS_UT_est_jobs.sh > CUR_submit.sh
    qsub CUR_submit.sh
done

rm -r "$TMPDIR"/VARsim_UT/
rm -r CUR_submit.sh



