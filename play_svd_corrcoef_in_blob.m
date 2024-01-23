clear all
close all
clc

root = 'C:\Users\User\Documents\fMRI_EXP\Alon\';
nVoxels = 117;
maskName = ['EC_HexOnHexPeak_mask_' num2str(nVoxels) 'vox'];
addpath(fullfile(root,'from_git','util'))
load(fullfile(root,'subspaceGener',['AUC_visualisation_' maskName '.mat']))

% subjects: sub-01 to sub-28
subjects = cell(1,28);
for iSub = 1:length(subjects)
    subjects{iSub}= ['sub-' num2str(iSub,'%02.f')];
end

% which betas: (for later analysis)
nPiles= 10;
nMaps = 4;
nRuns = 4;

nPC = 10;
cumVar_allSubs = nan(nMaps,nMaps,nRuns,nVoxels,length(subjects)); % nMapToBeProjected (single run) x nMap to get EVs from (3 runs of) x ...
H_all = nan(16,length(subjects));
for iSub = 1: length(subjects)

    betasToUse = squeeze(betasAllSubj(:,:,:,iSub));
    nVoxels = length(betasToUse(1,:,:));
    betaMeanLOO = zeros(40,4,nVoxels);

    nr = 1:4;

    %% LOO averages
    for r = 1:4
        betaMeanLOO(:,r,:) = mean(betasToUse(:,nr~=r,:),2); % mean of the betas in 3 of the 4 runs
    end

    Cmat = zeros(nMaps * nMaps, 4);

    for r = 1:4

        % Do SVD on the data averaged from the 3 runs except current run.
        % U is the spatial components, S the scalings (eigenvalues), V the
        % temporal components.

        %hex 1
        hex1avLoo = squeeze(betaMeanLOO(1:10,r,:));
        hex1avLoo = hex1avLoo - repmat(mean(hex1avLoo),nPiles,1); % substract the mean over condition
        [U1,S1,~] = svd(corr(hex1avLoo, hex1avLoo),'econ');%
        Cor_hex1_clean = U1(:, 1:nPC) * S1(1:nPC, 1:nPC) * U1(:, 1:nPC)';

        % hex 2
        hex2avLoo = squeeze(betaMeanLOO(11:20,r,:));
        hex2avLoo = hex2avLoo - repmat(mean(hex2avLoo),nPiles,1);
        [U2,S2,~] = svd(corr(hex2avLoo, hex2avLoo),'econ');%
        Cor_hex2_clean = U2(:, 1:nPC) * S2(1:nPC, 1:nPC) * U2(:, 1:nPC)';

        % cluster 1:
        clus1avLoo = squeeze(betaMeanLOO(21:30,r,:));
        clus1avLoo = clus1avLoo - repmat(mean(clus1avLoo),nPiles,1);
        [U3,S3,~] = svd(corr(clus1avLoo, clus1avLoo),'econ');%
        Cor_cluster1_clean = U3(:, 1:nPC) * S3(1:nPC, 1:nPC) * U3(:, 1:nPC)';

        % cluster 2:
        clus2avLoo = squeeze(betaMeanLOO(31:40,r,:));
        clus2avLoo = clus2avLoo - repmat(mean(clus2avLoo),nPiles,1);
        [U4,S4,~] = svd(corr(clus2avLoo, clus2avLoo),'econ');%
        Cor_cluster2_clean = U4(:, 1:nPC) * S4(1:nPC, 1:nPC) * U4(:, 1:nPC)';

        % Calculate SVD on the data from the run itself.
        % Sr are the scalings (eigenvalues)

        %hex 1
        hex1run = squeeze(betasToUse(1:10,r,:));
        hex1run =  hex1run - repmat(mean(hex1run,1),nPiles,1);
        [U1r,S1r,~] = svd(corr(hex1run, hex1run),'econ');
        Cor_hex1run_clean = U1r(:, 1:nPC) * S1r(1:nPC, 1:nPC) * U1r(:, 1:nPC)';

        %hex 2
        hex2run = squeeze(betasToUse(11:20,r,:));
        hex2run = hex2run - repmat(mean(hex2run,1),nPiles,1);
        [U2r,S2r,~] = svd(corr(hex2run, hex2run),'econ');
        Cor_hex2run_clean = U2r(:, 1:nPC) * S2r(1:nPC, 1:nPC) * U2r(:, 1:nPC)';

        %cluster 1:
        clus1run = squeeze(betasToUse(21:30,r,:));
        clus1run = clus1run - repmat(mean(clus1run),nPiles,1);
        [U3r,S3r,~] = svd(corr(clus1run, clus1run),'econ');
        Cor_cluster1run_clean = U3r(:, 1:nPC) * S3r(1:nPC, 1:nPC) * U3r(:, 1:nPC)';

        %cluster 2:
        clus2run = squeeze(betasToUse(31:40,r,:));
        clus2run = clus2run - repmat(mean(clus2run,1),nPiles,1);
        [U4r,S4r,~] = svd(corr(clus2run, clus2run),'econ');
        Cor_cluster2run_clean = U4r(:, 1:nPC) * S4r(1:nPC, 1:nPC) * U4r(:, 1:nPC)';

        %correlation of correlations matrices:
        % On hex 1:
        Cmat(1, r) = corr(Cor_hex1_clean(:), Cor_hex1run_clean(:));
        Cmat(2, r) = corr(Cor_hex1_clean(:), Cor_hex2run_clean(:));
        Cmat(3, r) = corr(Cor_hex1_clean(:), Cor_cluster1run_clean(:));
        Cmat(4, r) = corr(Cor_hex1_clean(:), Cor_cluster2run_clean(:));

        % On hex 2:
        Cmat(5, r) = corr(Cor_hex2_clean(:), Cor_hex1run_clean(:));
        Cmat(6, r) = corr(Cor_hex2_clean(:), Cor_hex2run_clean(:));
        Cmat(7, r) = corr(Cor_hex2_clean(:), Cor_cluster1run_clean(:));
        Cmat(8, r) = corr(Cor_hex2_clean(:), Cor_cluster2run_clean(:));

        % To cluster 1:
        Cmat(9, r) = corr(Cor_cluster1_clean(:), Cor_hex1run_clean(:));
        Cmat(10, r) = corr(Cor_cluster1_clean(:), Cor_hex2run_clean(:));
        Cmat(11, r) = corr(Cor_cluster1_clean(:), Cor_cluster1run_clean(:));
        Cmat(12, r) = corr(Cor_cluster1_clean(:), Cor_cluster2run_clean(:));

        % To cluster 2:
        Cmat(13, r) = corr(Cor_cluster2_clean(:), Cor_hex1run_clean(:));
        Cmat(14, r) = corr(Cor_cluster2_clean(:), Cor_hex2run_clean(:));
        Cmat(15, r) = corr(Cor_cluster2_clean(:), Cor_cluster1run_clean(:));
        Cmat(16, r) = corr(Cor_cluster2_clean(:), Cor_cluster2run_clean(:));

    end

    % get the average of the correlation

    H = mean(atanh(Cmat), 2);

    H_all(:, iSub) = H;
end

Hmean = mean(H_all, 2);
Hmat = reshape(Hmean, 4, 4);
imagesc(Hmat)
colorbar
