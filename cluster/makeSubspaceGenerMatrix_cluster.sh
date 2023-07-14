#!/bin/bash

sub=$1
scriptsDir=/home/fs0/abaram/scratch/shirley/alon/code/cluster
cat ${scriptsDir}/makeSubspaceGenerMatrix_cluster.m | sed "s/sub-00/${sub}/g" > ${scriptsDir}/${sub}_sl.m;
jobID1=`fsl_sub -N ${sub} -q long.q matlab -nodisplay -nosplash \< ${scriptsDir}/${sub}_sl.m`;
echo $jobID1