root = '/home/fs0/abaram/scratch/shirley/alon';

addpath(genpath(fullfile(root,'code')))

sub = '';
glm = 'glm1';
slName = 'L100';

makeSubspaceGenerMatrix_cor_cor(root,sub,glm,slName)
normAndSmoothSubspaceGenerMat_cor_cor(root,sub,slName)