function H = calcPileSimilarity(rawData,SPM)
% Shirley called this function pile_similarity_RSAtoolbox

% rawData: raw data from searchlight

% pre - whitening and calculate beta whights:
% default: normalize using all runs (overall):
% if wanted to be changed should add in the options: 'normmode':   
% 'overall': Does the multivariate noise normalisation overall (default)
% 'runwise': Does the multivariate noise normalisation by run

% take out Nan or inf:
mRowData = mean(rawData);
isNanOrInf = isnan(mRowData) + isinf(mRowData);
RowData2 = rawData(:,isNanOrInf==0);
[a,b]=size(rawData);
[a2,b2]=size(rawData);
if a~=a2||b~=b2
    disp('in calcPileSimilarity, check inf and nan, with/without data size: \n');
    disp(size(rawData));
    disp(size(RowData2));
end

[betaAll00,resMS,Sw_hat,beta_hat,shrinkage,trRR]=rsa.spm.noiseNormalizeBeta(RowData2,SPM);%,varargin);

% % structure of distance matrix:
DistMs = createDistStNew();

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

c1L = zeros(3,4);
c2L = zeros(3,4);
c3L = zeros(3,4);
c4L = zeros(3,4);

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
    rdm_hex1(r, :, :) = 1 - CC1;
    CC1 = atanh(CC1);

    c1L(1,r) = mean(CC1(DistMs.Hex==1)) / std(CC1(DistMs.Hex==1));
    c1L(2,r) = mean(CC1(DistMs.Hex==2)) / std(CC1(DistMs.Hex==2));
    c1L(3,r) = mean(CC1(DistMs.Hex==3)) / std(CC1(DistMs.Hex==3)); 
    
   % To hex 2:
    CC2  = corr(conSHex2av', conSHex2r');
    rdm_hex2(r, :, :) = 1 - CC2;
    CC2 = atanh(CC2);

    c2L(1,r) = mean(CC2(DistMs.Hex==1)) / std(CC2(DistMs.Hex==1));
    c2L(2,r) = mean(CC2(DistMs.Hex==2)) / std(CC2(DistMs.Hex==2));
    c2L(3,r) = mean(CC2(DistMs.Hex==3)) / std(CC2(DistMs.Hex==3));
      
    % cluster 1:
    CC3  = corr(conSC1av', conSC1r');
    rdm_cluster1(r, :, :) = CC3;
    CC3 = atanh(CC3);

    c3L(1,r) = mean(CC3(DistMs.ClustSmall==1)) / std(CC3(DistMs.ClustSmall==1));
    c3L(2,r) = mean(CC3(DistMs.ClustSmall==2)) / std(CC3(DistMs.ClustSmall==2)) ;
    c3L(3,r) = mean(CC3(DistMs.ClustSmall==3)) / std(CC3(DistMs.ClustSmall==3));
    
    % cluster 2
    CC4  = corr(conSC2av', conSC2r');
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

