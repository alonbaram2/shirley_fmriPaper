root = '/home/fs0/abaram/scratch/shirley/alon';

addpath(genpath(fullfile(root,'code')))

sub = '';
glm = 'glm1';
slName = 'L100';

makePileSimilarityRDM(root,sub,glm,slName)
normAndSmoothPileSimilarityRDM(root,sub,slName)