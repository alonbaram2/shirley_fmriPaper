root = '/home/fs0/abaram/scratch/shirley/alon';

addpath(genpath(fullfile(root,'code')))

sub = 'sub-27';
glm = 'glm1';
slName = 'L100';

makePileSimilarityRDM(root,sub,glm,slName)
normAndSmoothPileSimilarityRDM(root,sub,slName)