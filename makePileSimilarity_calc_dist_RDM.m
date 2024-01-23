function makePileSimilarity_calc_dist_RDM(root,sub,glm,slName)

glmDir = fullfile(root,'glms',glm,sub);
preprocDir = fullfile(root,'preproc',sub);
statsDir = fullfile(root,'RSA_pileSim_calculated_dist',sub);
mkdir(statsDir);

spm_path = '/home/fs0/abaram/scratch/MATLAB/spm12';
addpath(genpath(spm_path));
util_path = fullfile(root,'code','util');
addpath(genpath(util_path));

runs={'run-01','run-02','run-03','run-04'}; % run names (AKA sessions)
nRuns = length(runs);

% number of output files - 4x4 matrix where element i,j is the projection
% of data from condition i on eigenvalues fitted to condition j
nFiles = 8;

% searchlight file name
slDir = fullfile(root,'sl',sub);

% directory for temporary files:
tempFilesDir = fullfile(statsDir,'tempFiles');
mkdir(tempFilesDir);

% load SPM.mat:
load(fullfile(glmDir,'SPM.mat'));

% the names of the images created by the searchlight:  
outFiles = cell(1,nFiles);
for f=1:nFiles
    outFiles{f} = [statsDir,'/',slName,'_pileSim_calc_dist',num2str(f),'.nii']; 
end

% loading images:
inFiles = [];
for iRun = 1:nRuns
    run = runs{iRun};
    allfiles = cellstr(spm_select('FPList', [preprocDir,'/',run,'/'], '^r_cleaned_smoothed_.*.nii$')); % pre warping!
    inFiles = [inFiles; allfiles];
end

% loading searchlight definitions:
load(fullfile(slDir,[slName,'.mat']));

% run analysis. the details of the computation are in 
runSearchlight_alon(L,inFiles,outFiles,@pile_similarity_calc_dist_RSAtoolbox,run,tempFilesDir,'optionalParams',{SPM});