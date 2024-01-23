clear all
close all
clc

root = 'C:\Users\User\Documents\fMRI_EXP\Alon\';
nVoxels = 117;
maskName = ['EC_HexOnHexPeak_mask_' num2str(nVoxels) 'vox'];
addpath(fullfile(root,'from_git','util'))
load(fullfile(root,'subspaceGener',['AUC_visualisation_' maskName '.mat']))

% structure of distance matrix:
DistMs = creatDistStNew();


% subjects: sub-01 to sub-28
subjects = cell(1,28);
for iSub = 1:length(subjects)
    subjects{iSub}= ['sub-' num2str(iSub,'%02.f')];
end

% which betas: (for later analysis)
nPiles= 10;
nMaps = 4;
nRuns = 4;

nPC = 7;
lvvoxels = nVoxels;
cumVar_allSubs = nan(nMaps,nMaps,nRuns,nVoxels,length(subjects)); % nMapToBeProjected (single run) x nMap to get EVs from (3 runs of) x ...
H_all = nan(16,length(subjects));
for iSub = 1: length(subjects)

    betasToUse = squeeze(betasAllSubj(:,:,:,iSub));
    betaMeanLOO = zeros(40,4,nVoxels);

    nr = 1:4;

    %% LOO averages
    for r = 1:4
        betaMeanLOO(:,r,:) = mean(betasToUse(:,nr~=r,:),2); % mean of the betas in 3 of the 4 runs
    end
    rdm_hex1 = zeros(4, 10, 10);
    rdm_hex2 = zeros(4, 10, 10);
    rdm_cluster1 = zeros(4, 10, 10);
    rdm_cluster2 = zeros(4, 10, 10);

    % nPile = 10;
    for r = 1:4

        %averages:
        %hex 1
        hex1avLoo = squeeze(betaMeanLOO(1:10,r,:));
        hex1avLoo = hex1avLoo - repmat(mean(hex1avLoo),nPiles,1); % substract the mean over condition
        [U1,S1,~] = svd(hex1avLoo'*hex1avLoo,'econ');%
        hex1avLoo_clean = U1(:, 1:nPC)' * hex1avLoo';
        conSHex1av = hex1avLoo_clean - repmat(mean(hex1avLoo_clean),nPC,1);

        % hex 2
        hex2avLoo  = squeeze(betaMeanLOO(11:20,r,:));
        hex2avLoo = hex2avLoo - repmat(mean(hex2avLoo),nPiles,1); % substract the mean over condition
        [U2,S2,~] = svd(hex2avLoo'*hex2avLoo,'econ');%
        hex2avLoo_clean = U2(:, 1:nPC)' * hex2avLoo';
        conSHex2av =  hex2avLoo_clean - repmat(mean(hex2avLoo_clean), nPC, 1);

        % cluster 1:
        clus1avLoo = squeeze(betaMeanLOO(21:30,r,:));
        clus1avLoo = clus1avLoo - repmat(mean(clus1avLoo),nPiles,1);
        [U3,~,~] = svd(clus1avLoo'*clus1avLoo,'econ');%
        clus1avLoo_clean = U3(:, 1:nPC)' * clus1avLoo';
        conSC1av = clus1avLoo_clean - repmat(mean(clus1avLoo_clean), nPC, 1);

        % cluster 2:
        clus2avLoo = squeeze(betaMeanLOO(31:40,r,:));
        clus2avLoo = clus2avLoo - repmat(mean(clus2avLoo),nPiles,1);
        [U4,~,~] = svd(clus2avLoo'*clus2avLoo,'econ');%
        clus2avLoo_clean = U4(:, 1:nPC)' * clus2avLoo';
        conSC2av = clus2avLoo_clean - repmat(mean(clus2avLoo_clean), nPC, 1);

        % run patterns:
        %hex 1
        hex1run = squeeze(betasToUse(1:10, r,:));
        hex1run =  hex1run - repmat(mean(hex1run,1),nPiles,1);
        [U1r,S1r,~] = svd(hex1run'*hex1run,'econ');
        hex1run_clean = U1r(:, 1:nPC)' * hex1run';
        conSHex1r =  hex1run_clean - repmat(mean(hex1run_clean),nPC,1);

        %hex 2
        hex2run = squeeze(betasToUse(11:20, r,:));
        hex2run =  hex2run - repmat(mean(hex2run,1),nPiles,1);
        [U2r,S2r,~] = svd(hex2run'*hex2run,'econ');
        hex2run_clean = U2r(:, 1:nPC)' * hex2run';
        conSHex2r =  hex2run_clean - repmat(mean(hex2run_clean), nPC, 1);

        %cluster 1:
        clus1run = squeeze(betasToUse(21:30, r,:));
        clus1run =  clus1run - repmat(mean(clus1run,1),nPiles,1);
        [U3r,S3r,~] = svd(clus1run'*clus1run,'econ');
        clus1run_clean = U3r(:, 1:nPC)' * clus1run';
        conSC1r = clus1run_clean - repmat(mean(clus1run_clean), nPC, 1);

        %cluster 2:
        clus2run = squeeze(betasToUse(31:40, r,:));
        clus2run =  clus2run - repmat(mean(clus2run,1),nPiles,1);
        [U4r,S4r,~] = svd(clus2run'*clus2run,'econ');
        clus2run_clean = U4r(:, 1:nPC)' * clus2run';
        conSC2r = clus2run_clean - repmat(mean(clus2run_clean),nPC, 1);

        %correlations:

        % hex 1:
        CC1  = corr(conSHex1av, conSHex1r);
        rdm_hex1(r, :, :) = 1 - CC1;
        CC1 = atanh(CC1);

        c1L(1,r) = mean(CC1(DistMs.Hex==1)) / std(CC1(DistMs.Hex==1));
        c1L(2,r) = mean(CC1(DistMs.Hex==2)) / std(CC1(DistMs.Hex==2));
        c1L(3,r) = mean(CC1(DistMs.Hex==3)) / std(CC1(DistMs.Hex==3));

        % To hex 2:
        CC2  = corr(conSHex2av, conSHex2r);
        rdm_hex2(r, :, :) = 1 - CC2;
        CC2 = atanh(CC2);

        c2L(1,r) = mean(CC2(DistMs.Hex==1)) / std(CC2(DistMs.Hex==1));
        c2L(2,r) = mean(CC2(DistMs.Hex==2)) / std(CC2(DistMs.Hex==2));
        c2L(3,r) = mean(CC2(DistMs.Hex==3)) / std(CC2(DistMs.Hex==3));

        % cluster 1:
        CC3  = corr(conSC1av, conSC1r);
        rdm_cluster1(r, :, :) = CC3;
        CC3 = atanh(CC3);

        c3L(1,r) = mean(CC3(DistMs.ClustSmall==1)) / std(CC3(DistMs.ClustSmall==1));
        c3L(2,r) = mean(CC3(DistMs.ClustSmall==2)) / std(CC3(DistMs.ClustSmall==2)) ;
        c3L(3,r) = mean(CC3(DistMs.ClustSmall==3)) / std(CC3(DistMs.ClustSmall==3));

        % cluster 2
        CC4  = corr(conSC2av, conSC2r);
        rdm_cluster2(r, :, :) = CC4;
        CC4 = atanh(CC4);

        c4L(1,r) = mean(CC4(DistMs.ClustBig==1)) / std(CC4(DistMs.ClustBig==1));
        c4L(2,r) = mean(CC4(DistMs.ClustBig==2)) / std(CC4(DistMs.ClustBig==2)) ;
        c4L(3,r) = mean(CC4(DistMs.ClustBig==3)) / std(CC4(DistMs.ClustBig==3));


    end
    ms1 = mean(c1L,2);
    ms2 = mean(c2L,2);
    ms3 = mean(c3L,2);
    ms4 = mean(c4L,2);

    mean_rdm_hex1 = mean(rdm_hex1);
    mean_rdm_hex2 = mean(rdm_hex2);
    mean_rdm_cluster1 = mean(rdm_cluster1);
    mean_rdm_cluster2 = mean(rdm_cluster2);

    kandell_similarity_hex1 = corr(mean_rdm_hex1(DistMs.Hex > 0), DistMs.Hex(DistMs.Hex > 0) ,'type','Kendall');
    kandell_similarity_hex2 = corr(mean_rdm_hex2(DistMs.Hex > 0), DistMs.Hex(DistMs.Hex > 0) ,'type','Kendall');
    kandell_similarity_cluster1 = corr(mean_rdm_cluster1(DistMs.ClustSmall > 0), DistMs.ClustSmall(DistMs.ClustSmall > 0),'type','Kendall');
    kandell_similarity_cluster2 = corr(mean_rdm_cluster2(DistMs.ClustBig > 0), DistMs.ClustBig(DistMs.ClustBig > 0),'type','Kendall');

    H = [ms1;ms2;ms3;ms4;kandell_similarity_hex1;kandell_similarity_hex2;kandell_similarity_cluster1;kandell_similarity_cluster2];


    H_all(:, iSub) = H;

end

H_all_mean = mean(H_all, 2);