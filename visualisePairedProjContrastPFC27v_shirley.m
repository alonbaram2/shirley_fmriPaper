clear
close all
clc

exp_root = 'C:\Users\User\Documents\fMRI_EXP\';
root = 'C:\Users\User\Documents\fMRI_EXP\Alon\';
dataDir = fullfile(root,"beh","training");
%addpath(genpath(fullfile(root,'code')));
addpath(genpath("C:\Users\User\OneDrive\Desktop\spm12"))

subjects = cell(1,28);
for iSub = 1:length(subjects)
    subjects{iSub}= ['sub-' num2str(iSub,'%02.f')];   
end


blob_corrdinates = [-4,44,-20];% Alon PFC blob coordinate

one_size = 3;
titleStr = ['PFC', num2str(2*one_size + 1)]
allElem = 1:16;
proHexElem = [1, 2, 5, 6];
proHexClElem = [9, 10, 13, 14];
proClClElem = [11, 12, 15, 16];
proClHexElem = [3, 4, 7, 8];
allData = nan(length(allElem),length(subjects));
for iSub=1:length(subjects)
    sub = subjects{iSub};
    subspaceGenerDir = fullfile(root,'subspaceGener',sub);

    % Full matrix
    for iElem = 1:16
        fname{iElem} = fullfile(subspaceGenerDir,['smth6mni_L100_projMat' num2str(allElem(iElem)) '.nii']);
        V_all = spm_vol(fname{iElem});
        [ROI_dat1,XYZ1] = spm_read_vols(V_all,0);
        x = XYZ1(1,:);
        y = XYZ1(2,:);
        z = XYZ1(3,:);
        ind27 = zeros(size(x));% indexes for voxels in the blob
        for jx = -one_size :one_size %-4:2:4
            for jy = -one_size :one_size %-4:2:4
                for jz = -one_size :one_size % -4:2:4
                    ind27 = ind27 +(x==(blob_corrdinates(1)+jx)).*(y==(blob_corrdinates(2)+jy)).*(z==(blob_corrdinates(3)+jz));
                end
            end
        end
        allData(iElem,iSub) = mean(ROI_dat1(ind27==1));          
    end
end

proHex_allData = allData(proHexElem,:);
proHexCl_allData = allData(proHexClElem,:);
proClCl_allData = allData(proClClElem, :);
proClHex_allData = allData(proClHexElem, :);

mean_matrixHex_suqres = mean(proHex_allData);
mean_matrixHexCl_suqres = mean(proHexCl_allData);
mean_matrixClCl_suqres = mean(proClCl_allData);
mean_matrixClHex_suqres = mean(proClHex_allData);

[h_hex_cl, p_hex_cl] = ttest(mean_matrixHex_suqres, mean_matrixHexCl_suqres, 'Tail', 'right');
[h_cl_hex, p_cl_hex] = ttest(mean_matrixHex_suqres, mean_matrixClHex_suqres, 'Tail', 'right');
[h_clcl_hex, p_clcl_hex, ~, stats_hex] = ttest(mean_matrixClCl_suqres, mean_matrixHexCl_suqres, 'Tail', 'right');
[h_clcl_clhex, p_clcl_clhex, ~, stats] = ttest(mean_matrixClCl_suqres, mean_matrixClHex_suqres, 'Tail', 'right');
[h_clcl_clhex_sym, p_clcl_clhex_sym] = ttest(mean_matrixClCl_suqres, mean_matrixClHex_suqres);

cl_contrast = mean_matrixClCl_suqres - mean_matrixClHex_suqres;

load(fullfile(dataDir,'pConCsCcor2.mat'))
load(fullfile(dataDir,'pConCsCncor2.mat'))
load(fullfile(dataDir,'pConCsCncor1.mat'))

cor_minus_notcor = pConCsCcor2 - pConCsCncor2;

[c, p] = corrcoef(cl_contrast, cor_minus_notcor);
[c1, p1] = corrcoef(cl_contrast, pConCsCncor1);
[c2, p2] = corrcoef(cl_contrast, pConCsCncor2);

projMatAllSubj = reshape(allData,[4,4,length(subjects)]);

projMatAllSubj2 = projMatAllSubj;
projMatAllSubj2(: ,4 ,:) = projMatAllSubj(:, 3, :);
projMatAllSubj2(: ,3 ,:) = projMatAllSubj(:, 4, :);
v = projMatAllSubj2(4, :, :);
projMatAllSubj2(4 ,: ,:) = projMatAllSubj2(3, :, :);
projMatAllSubj2(3 ,: ,:) = v;

figure
hold on
title([titleStr ', mean'])
projMatMean = squeeze(mean(projMatAllSubj2,3));
imagesc(projMatMean);
colormap default
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij

figure
hold on
title([titleStr ', t-value'])
[~,~,~,S] = ttest(projMatAllSubj2,50,'dim',3);
projMatT = S.tstat;
imagesc(projMatT);
colormap default
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij

figure(10)
subplot(1,2,1)
projMatMean = squeeze(mean(projMatAllSubj2,3));
imagesc(projMatMean);
title([titleStr ', mean'])
colormap default
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij

subplot(1,2,2)
[~,~,~,S] = ttest(projMatAllSubj2,50,'dim',3);
projMatT = S.tstat;
imagesc(projMatT);
title([titleStr ', t-value'])
colormap default
labels = {'Hl','Hs','Cl','Cs'};
xticklabels(labels);
yticklabels(labels);
xticks(1:4)
yticks(1:4)
ylabel('Data to project')
xlabel('Data to calculate subspace')
colorbar
axis ij