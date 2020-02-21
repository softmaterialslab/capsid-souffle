#!/bin/bash
for INPUT_TEMP in 149
do
dir="/N/dc2/scratch/lm44/MVM_T_$INPUT_TEMP"
mkdir $dir
mkdir $dir/infiles
mkdir $dir/outfiles
cp -r bin/infiles/ $dir
cp bin/Makefile $dir
cp src/capsid-souffle $dir
cp scripts/iu_cluster_job_script.pbs $dir
cd $dir
sed -i 's/INPUT_TEMP/'$INPUT_TEMP'/g' iu_cluster_job_script.pbs
sbatch iu_cluster_job_script.pbs
cd
cd ~/capsid-souffle-fork
done
