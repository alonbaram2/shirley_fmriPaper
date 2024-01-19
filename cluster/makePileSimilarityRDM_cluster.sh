#!/bin/bash

sub=$1
scriptsDir=/home/fs0/abaram/scratch/shirley/alon/code/cluster
cat ${scriptsDir}/makePileSimilarityRDM_cluster.m | sed "s/sub-00/${sub}/g" > ${scriptsDir}/${sub}_PileSimilarity.m;
jobID1=`fsl_sub -N ${sub}_pileSimilarity -q long.q matlab -nodisplay -nosplash \< ${scriptsDir}/${sub}_PileSimilarity.m`;
echo $jobID1