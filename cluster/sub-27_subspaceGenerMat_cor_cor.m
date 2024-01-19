root = '/home/fs0/abaram/scratch/shirley/alon';

addpath(genpath(fullfile(root,'code')))

sub = 'sub-27';
glm = 'glm1';
slName = 'L100';

makeSubspaceGenerMatrix_cor_cor(root,sub,glm,slName)
normAndSmoothSubspaceGenerMat_cor_cor(root,sub,slName)