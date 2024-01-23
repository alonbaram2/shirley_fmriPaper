function H = pile_similarity_calc_dist_RSAtoolbox(RowData,SPM)
% RowData: raw data from searchlight

% pre - whitening and calculate beta whights:
% default: normalize using all runs (overall):
% if wanted to be changed should add in the options: 'normmode':   
% 'overall': Does the multivariate noise normalisation overall (default)
% 'runwise': Does the multivariate noise normalisation by run

% take out Nan or inf:
mRowData = mean(RowData);
isNanOrInf = isnan(mRowData) + isinf(mRowData);
RowData2 = RowData(:,isNanOrInf==0);
[a,b]=size(RowData);
[a2,b2]=size(RowData);
if a~=a2||b~=b2
    disp('in pile_similarity_RSAtoolbox, check inf and nan, with/without data size: \n');
    disp(size(RowData));
    disp(size(RowData2));
end

[betaAll00,resMS,Sw_hat,beta_hat,shrinkage,trRR]=noiseNormalizeBeta(RowData2,SPM);%,varargin);

% % structure of distance matrix:
load('calc_distance_matrixes_struct.mat');

% which betas: (for later analysis)
numPile= 10;
nMap = 4;
v = ones(1,numPile*nMap);
zeroBetas = 10 + 5 +2 +7;
bn0 =  [zeros(1,5) v zeros(1,zeroBetas)];%zeros in the end - betas for 5 maps,the catch, the start and the movement regressors + ventricles regressors
bn = [bn0 bn0 bn0 bn0];
nwb = length(bn);
b_idx0 = 1:nwb;
b_idx = b_idx0(bn==1);

betaAll0 = betaAll00(b_idx,:);

lvvoxels = length(betaAll0(1,:));
betaAll1 = reshape(betaAll0,[40,4,lvvoxels]);
betaAllLooAv = zeros(40,4,lvvoxels);

nr = 1:4;
%% LOO averages:
for r = 1:4
    betaAllLooAv(:,r,:) = mean(betaAll1(:,nr~=r,:),2);
end

rdm_hex1 = zeros(4, 10, 10);
rdm_hex2 = zeros(4, 10, 10);
rdm_cluster1 = zeros(4, 10, 10);
rdm_cluster2 = zeros(4, 10, 10);

% nPile = 10;
for r = 1:4
    
    %averages:
    %hex 1
    conSHex1av = squeeze(betaAllLooAv(1:10,r,:));
    conSHex1av = conSHex1av - repmat(mean(conSHex1av,2),1,lvvoxels);
    
    % hex 2
    conSHex2av = squeeze(betaAllLooAv(11:20,r,:));
    conSHex2av = conSHex2av - repmat(mean(conSHex2av,2),1,lvvoxels);
    
    % cluster 1:
    conSC1av = squeeze(betaAllLooAv(21:30,r,:));
    conSC1av = conSC1av - repmat(mean(conSC1av,2),1,lvvoxels); 
    
    % cluster 2:
    conSC2av = squeeze(betaAllLooAv(31:40,r,:)); 
    conSC2av = conSC2av - repmat(mean(conSC2av,2),1,lvvoxels); 
    
    % run patterns:
    %hex 1
    conSHex1r = betaAll0((r-1)*40+1:(r-1)*40+10,:);
    conSHex1r =  conSHex1r - repmat(mean(conSHex1r,2),1,lvvoxels);
    
    %hex 2
    conSHex2r = betaAll0((r-1)*40+11:(r-1)*40+20,:);
    conSHex2r =  conSHex2r - repmat(mean(conSHex2r,2),1,lvvoxels);
    
    %cluster 1:
    conSC1r = betaAll0((r-1)*40+21:(r-1)*40+30,:);
    conSC1r = conSC1r - repmat(mean(conSC1r,2),1,lvvoxels);
    
    %cluster 2:
    conSC2r = betaAll0((r-1)*40+31:(r-1)*40+40,:);
    conSC2r = conSC2r - repmat(mean(conSC2r,2),1,lvvoxels);
    
    %correlations:

    % hex 1:    
    CC1  = corr(conSHex1av', conSHex1r');
    rdm_hex_big(r, :, :) = 1 - CC1;
    
   % To hex 2:
    CC2  = corr(conSHex2av', conSHex2r');
    rdm_hex_small(r, :, :) = 1 - CC2;

    % cluster 1:
    CC3  = corr(conSC1av', conSC1r');
    rdm_cluster_small(r, :, :) = CC3;
    
    % cluster 2
    CC4  = corr(conSC2av', conSC2r');
    rdm_cluster_big(r, :, :) = CC4;

end

mean_rdm_hex_big = mean(rdm_hex_big);
mean_rdm_hex_small = mean(rdm_hex_small);
mean_rdm_cluster_small = mean(rdm_cluster_small);
mean_rdm_cluster_big = mean(rdm_cluster_big);

% correct model rdm:
kandell_similarity_hex_big = corr(mean_rdm_hex1(distance_matrixes_struct.pile_dist_hex.Hex > 0), distance_matrixes_struct.pile_dist_hex(distance_matrixes_struct.pile_dist_hex > 0) ,'type','Kendall');
kandell_similarity_hex_small = corr(mean_rdm_hex2(distance_matrixes_struct.pile_dist_hex.Hex > 0), distance_matrixes_struct.pile_dist_hex(distance_matrixes_struct.pile_dist_hex > 0) ,'type','Kendall');
kandell_similarity_cluster_small = corr(mean_rdm_cluster_small(distance_matrixes_struct.pile_dist_cluster_small > 0), distance_matrixes_struct.pile_dist_cluster_small(istance_matrixes_struct.pile_dist_cluster_smalll > 0),'type','Kendall');
kandell_similarity_cluster_big = corr(mean_rdm_clusterbig(distance_matrixes_struct.pile_dist_cluster_big > 0), distance_matrixes_struct.pile_dist_cluster_big(distance_matrixes_struct.pile_dist_cluster_big > 0),'type','Kendall');

% wrong model rdm:
kandell_similarity_hex_cluster_big = corr(mean_rdm_hex1(distance_matrixes_struct.pile_dist_hex.Hex > 0),  distance_matrixes_struct.pile_dist_cluster_big(distance_matrixes_struct.pile_dist_cluster_big > 0) ,'type','Kendall');
kandell_similarity_hex_cluster_small = corr(mean_rdm_hex2(distance_matrixes_struct.pile_dist_hex.Hex > 0), distance_matrixes_struct.pile_dist_cluster_small(istance_matrixes_struct.pile_dist_cluster_smalll > 0) ,'type','Kendall');
kandell_similarity_cluster_hex_small = corr(mean_rdm_cluster_small(distance_matrixes_struct.pile_dist_cluster_small > 0),distance_matrixes_struct.pile_dist_hex(distance_matrixes_struct.pile_dist_hex > 0),'type','Kendall');
kandell_similarity_cluster_hex_big = corr(mean_rdm_clusterbig(distance_matrixes_struct.pile_dist_cluster_big > 0), distance_matrixes_struct.pile_dist_hex(distance_matrixes_struct.pile_dist_hex > 0),'type','Kendall');

H = [kandell_similarity_hex_big;kandell_similarity_hex_small;kandell_similarity_cluster_small;kandell_similarity_cluster_big;kandell_similarity_hex_cluster_big; kandell_similarity_hex_cluster_small; kandell_similarity_cluster_hex_small; kandell_similarity_cluster_hex_big];

