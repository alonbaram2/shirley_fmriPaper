function H = calcSubspaceGener_cor_cor(rawData,SPM)
% rawData: raw data from searchlight
% calculate the subspace generalisation (instead of area under the curve of explained
% variance as a function of number of EVs calculate correlation of correlation matrices)

flagCal = 1; % flag whether to caclulate subspace generalisatoin, will change to 0 only if there is a problem with NaNs or Infs. 

% if there are any voxels in the SL that have data with NaN of Inf, get rid
% of them 
isNanOrInf = isnan(mean(rawData,1)) + isinf(mean(rawData,1));
rawData = rawData(:,isNanOrInf==0);

% pre - whitening and calculate betas:
% default: normalize using all runs (overall):
[betaAll00,~,~,~,~,~]=rsa.spm.noiseNormalizeBeta(rawData,SPM);

% which betas: (for later analysis)
nPiles= 10;
nMaps = 4;

pilesBetasInds = ones(1,nPiles*nMaps);
zeroBetas = 10 + 5 +2 +7; % 10: piles in map 5; 5: catch trials regs, one per block; 2: catch trial question and start experiment reg; 7: nuisance regressor (6 motion + CSF)
betasToUse_singleRun =  [zeros(1,5) pilesBetasInds zeros(1,zeroBetas)]; % the first 5 0s: The regs modelling the random walks at the start of each map. 
betasToUse = [betasToUse_singleRun betasToUse_singleRun betasToUse_singleRun betasToUse_singleRun];
% get indeces of betas to use
b_idx0 = 1:length(betasToUse);
b_idx = b_idx0(betasToUse==1);

% check there are voxels with no betas that are NaNs or Infs
betasToUse = betaAll00(b_idx,:);
isNanOrInf = isnan(mean(betasToUse,1))  + isinf(mean(betasToUse,1));

if sum(isNanOrInf>0)
    disp(' \n \n \n!!!!! There was a NaN or Inf, did not calculate !!!!! \n \n \n ')
    warning(' \n \n \n!!!!! There was a NaN or Inf, did not calculate !!!!! \n \n \n ')
    flagCal = 0;
end


if flagCal==1 % if betas are OK calculate, otherwise put -10000 at that voxel
    nVoxels = length(betasToUse(1,:));
    betasToUse = reshape(betasToUse,[40,4,nVoxels]); % nPiles*nMaps, nRuns, nVoxels
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
        Cor_hex1_clean = corr(hex1avLoo, hex1avLoo);
        
        % hex 2
        hex2avLoo = squeeze(betaMeanLOO(11:20,r,:));
        Cor_hex2_clean = corr(hex2avLoo, hex2avLoo);

        % cluster 1:
        clus1avLoo = squeeze(betaMeanLOO(21:30,r,:));
        Cor_cluster1_clean = corr(clus1avLoo, clus1avLoo);

        % cluster 2:
        clus2avLoo = squeeze(betaMeanLOO(31:40,r,:));
        Cor_cluster2_clean = corr(clus2avLoo, clus2avLoo);

        % Calculate SVD on the data from the run itself. 
        % Sr are the scalings (eigenvalues)
        
        %hex 1
        hex1run = squeeze(betasToUse(1:10,r,:));
        Cor_hex1run_clean = corr(hex1run, hex1run);

        %hex 2
        hex2run = squeeze(betasToUse(11:20,r,:));
        Cor_hex2run_clean = corr(hex2run, hex2run);

        %cluster 1:
        clus1run = squeeze(betasToUse(21:30,r,:));
        Cor_cluster1run_clean = corr(clus1run, clus1run);

        %cluster 2:
        clus2run = squeeze(betasToUse(31:40,r,:));
        Cor_cluster2run_clean = corr(clus2run, clus2run);

        %correlation of correlations matrices:
        % On hex 1:
        Cmat(1, r) = corr(atanh(Cor_hex1_clean(:)), atanh(Cor_hex1run_clean(:)));
        Cmat(2, r) = corr(atanh(Cor_hex1_clean(:)), atanh(Cor_hex2run_clean(:)));
        Cmat(3, r) = corr(atanh(Cor_hex1_clean(:)), atanh(Cor_cluster1run_clean(:)));
        Cmat(4, r) = corr(atanh(Cor_hex1_clean(:)), atanh(Cor_cluster2run_clean(:)));
        
        % On hex 2:
        Cmat(5, r) = corr(atanh(Cor_hex2_clean(:)), atanh(Cor_hex1run_clean(:)));
        Cmat(6, r) = corr(atanh(Cor_hex2_clean(:)), atanh(Cor_hex2run_clean(:)));
        Cmat(7, r) = corr(atanh(Cor_hex2_clean(:)), atanh(Cor_cluster1run_clean(:)));
        Cmat(8, r) = corr(atanh(Cor_hex2_clean(:)), atanh(Cor_cluster2run_clean(:)));
        
        % To cluster 1:
        Cmat(9, r) = corr(atanh(Cor_cluster1_clean(:)), atanh(Cor_hex1run_clean(:)));
        Cmat(10, r) = corr(atanh(Cor_cluster1_clean(:)), atanh(Cor_hex2run_clean(:)));
        Cmat(11, r) = corr(atanh(Cor_cluster1_clean(:)), atanh(Cor_cluster1run_clean(:)));
        Cmat(12, r) = corr(atanh(Cor_cluster1_clean(:)), atanh(Cor_cluster2run_clean(:)));
        
        % To cluster 2:
        Cmat(13, r) = corr(atanh(Cor_cluster2_clean(:)), atanh(Cor_hex1run_clean(:)));
        Cmat(14, r) = corr(atanh(Cor_cluster2_clean(:)), atanh(Cor_hex2run_clean(:)));
        Cmat(15, r) = corr(atanh(Cor_cluster2_clean(:)), atanh(Cor_cluster1run_clean(:)));
        Cmat(16, r) = corr(atanh(Cor_cluster2_clean(:)), atanh(Cor_cluster2run_clean(:)));
        
    end
    
    % get the average of the area under the curve across all runs

    H = mean(atanh(Cmat), 2);
    
else
    H = -10000*ones(1,16);
end
