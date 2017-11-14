#!/bin/bash

mkdir "$TMPDIR"/varSim/

cd "$HOME"/varSim

for i in `seq 1 100`;
do
	sed s/iter/$i/g VS_data_jobs.sh > CUR_submit.sh
    qsub CUR_submit.sh
done

rm -r "$TMPDIR"/varSim/
rm -r CUR_submit.sh



