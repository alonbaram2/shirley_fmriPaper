#!/bin/bash

sub=$1
scriptsDir=/home/fs0/abaram/scratch/shirley/alon/code/cluster
cat ${scriptsDir}/makeSubspaceGenerMatrix_cor_cor_cluster.m | sed "s/sub-00/${sub}/g" > ${scriptsDir}/${sub}_subspaceGenerMat_cor_cor.m;
jobID1=`fsl_sub -N ${sub}_Mat_cor_cor -q long.q matlab -nodisplay -nosplash \< ${scriptsDir}/${sub}_subspaceGenerMat_cor_cor.m`;
echo $jobID1