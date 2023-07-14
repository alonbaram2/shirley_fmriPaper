#!/bin/bash

root=/vols/Scratch/abaram/shirley/alon/code/cluster
for iSub in $(seq -f "%02g" 2 28)
do
  sub=sub-${iSub}
  echo $sub
  ${root}/makeSubspaceGenerMatrix_cluster.sh $sub
done