root = '/home/fs0/abaram/scratch/shirley/alon';

addpath(genpath(fullfile(root,'code')))

sub = 'sub-10';
glm = 'glm1';
slName = 'L100';

makeSubspaceGenerMatrix(root,sub,glm,slName)