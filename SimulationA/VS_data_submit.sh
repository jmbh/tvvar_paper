#!/bin/bash

mkdir "$TMPDIR"/tvvarSim2019/

cd "$HOME"/tvvarSim2019

for i in `seq 1 100`;
do
	sed s/iter/$i/g VS_data_jobs.sh > CUR_submit.sh
    sbatch CUR_submit.sh
done

rm -r "$TMPDIR"/tvvarSim2019/
rm -r CUR_submit.sh



