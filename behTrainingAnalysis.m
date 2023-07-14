close all
clear
clc

dataDir = '/Volumes/Scratch_abaram/shirley/alon/beh/training';
% dataDir = '/vols/Scratch/abaram/shirley/alon/beh/training';

load(fullfile(dataDir,'allTrainingData.mat'),'D');

nSub = length(D.middle.hex_day1_ratioCorrect);

figure

%% is in middle task

% plot all datapoints
subplot(2,3,1)
hold on
title('Is in middle?')
% Plot individual data points
scatter(ones(nSub,1),100*D.middle.hex_day1_ratioCorrect,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(2*ones(nSub,1),100*D.middle.hex_day2_ratioCorrect,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(3*ones(nSub,1),100*D.middle.clus_day3_ratioCorrect,3,'filled','b','MarkerFaceAlpha',0.5)
scatter(4*ones(nSub,1),100*D.middle.clus_day4_ratioCorrect,3,'filled','b','MarkerFaceAlpha',0.5)

% plot means and error bars
wBar = 0.15;
plot([1+0.2-wBar, 1+0.2+wBar], 100*[mean(D.middle.hex_day1_ratioCorrect),mean(D.middle.hex_day1_ratioCorrect)],'k','lineWidth',2)
plot([2-0.2-wBar, 2-0.2+wBar], 100*[mean(D.middle.hex_day2_ratioCorrect),mean(D.middle.hex_day2_ratioCorrect)],'k','lineWidth',2)
plot([3+0.2-wBar, 3+0.2+wBar], 100*[mean(D.middle.clus_day3_ratioCorrect),mean(D.middle.clus_day3_ratioCorrect)],'k','lineWidth',2)
plot([4-0.2-wBar, 4-0.2+wBar], 100*[mean(D.middle.clus_day4_ratioCorrect),mean(D.middle.clus_day4_ratioCorrect)],'k','lineWidth',2)

errorbar(1+0.2,100*mean(D.middle.hex_day1_ratioCorrect),100*std(D.middle.hex_day1_ratioCorrect)/sqrt(nSub),'k')
errorbar(2-0.2,100*mean(D.middle.hex_day2_ratioCorrect),100*std(D.middle.hex_day2_ratioCorrect)/sqrt(nSub),'k')
errorbar(3+0.2,100*mean(D.middle.clus_day3_ratioCorrect),100*std(D.middle.clus_day3_ratioCorrect)/sqrt(nSub),'k')
errorbar(4-0.2,100*mean(D.middle.clus_day4_ratioCorrect),100*std(D.middle.clus_day4_ratioCorrect)/sqrt(nSub),'k')

% Plot lines connecting the data points
for iSub=1:nSub
    plot([1,2],100*[D.middle.hex_day1_ratioCorrect(iSub),D.middle.hex_day2_ratioCorrect(iSub)],'color',[1 0 0 0.2])
    plot([3,4],100*[D.middle.clus_day3_ratioCorrect(iSub),D.middle.clus_day4_ratioCorrect(iSub)],'color',[0 0 1 0.2])
end
xlim([0,5])
labels = {'hex day1','hex day2','comm day3', 'comm day4'};
set(gca,'FontSize',14,'XTick',1:4,'XTickLabels',labels,'XTickLabelRotation', 45)
ylabel('% correct')


% stats for whether answers oneach day were different from chanve
[~, pMiddleHexDay1] = ttest(D.middle.hex_day1_ratioCorrect,0.5,'Tail','right')
[~, pMiddleHexDay2] = ttest(D.middle.hex_day2_ratioCorrect,0.5,'Tail','right')
[~, pMiddleClusDay3] = ttest(D.middle.clus_day3_ratioCorrect,0.5,'Tail','right')
[~, pMiddleClusDay4] = ttest(D.middle.clus_day4_ratioCorrect,0.5,'Tail','right')

% stats for improvement beween days
[~, pMiddleHexDays] = ttest2(D.middle.hex_day2_ratioCorrect,D.middle.hex_day1_ratioCorrect,'Tail','right')
[~, pMiddleClusDays] = ttest2(D.middle.clus_day4_ratioCorrect,D.middle.clus_day3_ratioCorrect,'Tail','right')

%% piles task
subplot(2,3,2)
hold on
title('Is this image next?')
% Plot individual data points
scatter(ones(nSub,1),100*D.piles.hex_day1_ratioCorrect,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(2*ones(nSub,1),100*D.piles.hex_day2_ratioCorrect,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(3*ones(nSub,1),100*D.piles.clus_day3_ratioCorrect,3,'filled','b','MarkerFaceAlpha',0.5)
scatter(4*ones(nSub,1),100*D.piles.clus_day4_ratioCorrect,3,'filled','b','MarkerFaceAlpha',0.5)

% plot means and SEM error bars
wBar = 0.15;
plot([1+0.2-wBar, 1+0.2+wBar], 100*[mean(D.piles.hex_day1_ratioCorrect),mean(D.piles.hex_day1_ratioCorrect)],'k','lineWidth',2)
plot([2-0.2-wBar, 2-0.2+wBar], 100*[mean(D.piles.hex_day2_ratioCorrect),mean(D.piles.hex_day2_ratioCorrect)],'k','lineWidth',2)
plot([3+0.2-wBar, 3+0.2+wBar], 100*[mean(D.piles.clus_day3_ratioCorrect),mean(D.piles.clus_day3_ratioCorrect)],'k','lineWidth',2)
plot([4-0.2-wBar, 4-0.2+wBar], 100*[mean(D.piles.clus_day4_ratioCorrect),mean(D.piles.clus_day4_ratioCorrect)],'k','lineWidth',2)

errorbar(1+0.2,100*mean(D.piles.hex_day1_ratioCorrect),100*std(D.piles.hex_day1_ratioCorrect)/sqrt(nSub),'k')
errorbar(2-0.2,100*mean(D.piles.hex_day2_ratioCorrect),100*std(D.piles.hex_day2_ratioCorrect)/sqrt(nSub),'k')
errorbar(3+0.2,100*mean(D.piles.clus_day3_ratioCorrect),100*std(D.piles.clus_day3_ratioCorrect)/sqrt(nSub),'k')
errorbar(4-0.2,100*mean(D.piles.clus_day4_ratioCorrect),100*std(D.piles.clus_day4_ratioCorrect)/sqrt(nSub),'k')

% Plot lines connecting the data points
for iSub=1:nSub
    plot([1,2],100*[D.piles.hex_day1_ratioCorrect(iSub),D.piles.hex_day2_ratioCorrect(iSub)],'color',[1 0 0 0.2])
    plot([3,4],100*[D.piles.clus_day3_ratioCorrect(iSub),D.piles.clus_day4_ratioCorrect(iSub)],'color',[0 0 1 0.2])
end

xlim([0,5]);
ylabel('% correct')
labels = {'hex day1','hex day2','comm day3', 'comm day4'};
set(gca,'FontSize',14,'XTick',1:4,'XTickLabels',labels,'XTickLabelRotation', 45)



% stats for whether answers oneach day were different from chanve
[~, pPilesHexDay1] = ttest(D.piles.hex_day1_ratioCorrect,0.5,'Tail','right')
[~, pPilesHexDay2] = ttest(D.piles.hex_day2_ratioCorrect,0.5,'Tail','right')
[~, pPilesClusDay3] = ttest(D.piles.clus_day3_ratioCorrect,0.5,'Tail','right')
[~, pPilesClusDay4] = ttest(D.piles.clus_day4_ratioCorrect,0.5,'Tail','right')

% stats for improvement beween days
[~, pPilesHexDays] = ttest2(D.piles.hex_day2_ratioCorrect,D.piles.hex_day1_ratioCorrect,'Tail','right')
[~, pPilesClusDays] = ttest2(D.piles.clus_day4_ratioCorrect,D.piles.clus_day3_ratioCorrect,'Tail','right')


%% which is closer task

subplot(2,3,3)
hold on
title('Which is closer?')
% Plot individual data points
scatter(ones(nSub,1),100*D.closer.hex_day1_ratioCorrect,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(2*ones(nSub,1),100*D.closer.hex_day2_ratioCorrect,3,'filled','r','MarkerFaceAlpha',0.5)
scatter(3*ones(nSub,1),100*D.closer.clus_day3_ratioCorrect,3,'filled','b','MarkerFaceAlpha',0.5)
scatter(4*ones(nSub,1),100*D.closer.clus_day4_ratioCorrect,3,'filled','b','MarkerFaceAlpha',0.5)

% plot means and SEM error bars
wBar = 0.15;
plot([1+0.2-wBar, 1+0.2+wBar], 100*[mean(D.closer.hex_day1_ratioCorrect),mean(D.closer.hex_day1_ratioCorrect)],'k','lineWidth',2)
plot([2-0.2-wBar, 2-0.2+wBar], 100*[mean(D.closer.hex_day2_ratioCorrect),mean(D.closer.hex_day2_ratioCorrect)],'k','lineWidth',2)
plot([3+0.2-wBar, 3+0.2+wBar], 100*[mean(D.closer.clus_day3_ratioCorrect),mean(D.closer.clus_day3_ratioCorrect)],'k','lineWidth',2)
plot([4-0.2-wBar, 4-0.2+wBar], 100*[mean(D.closer.clus_day4_ratioCorrect),mean(D.closer.clus_day4_ratioCorrect)],'k','lineWidth',2)

errorbar(1+0.2,100*mean(D.closer.hex_day1_ratioCorrect),100*std(D.closer.hex_day1_ratioCorrect)/sqrt(nSub),'k')
errorbar(2-0.2,100*mean(D.closer.hex_day2_ratioCorrect),100*std(D.closer.hex_day2_ratioCorrect)/sqrt(nSub),'k')
errorbar(3+0.2,100*mean(D.closer.clus_day3_ratioCorrect),100*std(D.closer.clus_day3_ratioCorrect)/sqrt(nSub),'k')
errorbar(4-0.2,100*mean(D.closer.clus_day4_ratioCorrect),100*std(D.closer.clus_day4_ratioCorrect)/sqrt(nSub),'k')

% Plot lines connecting the data points
for iSub=1:nSub
    plot([1,2],100*[D.closer.hex_day1_ratioCorrect(iSub),D.closer.hex_day2_ratioCorrect(iSub)],'color',[1 0 0 0.2])
    plot([3,4],100*[D.closer.clus_day3_ratioCorrect(iSub),D.closer.clus_day4_ratioCorrect(iSub)],'color',[0 0 1 0.2])
end

xlim([0,5]);
ylabel('% correct')
labels = {'hex day1','hex day2','comm day3', 'comm day4'};
set(gca,'FontSize',14,'XTick',1:4,'XTickLabels',labels,'XTickLabelRotation', 45)



% stats for whether answers oneach day were different from chanve
[~, pCloserHexDay1] = ttest(D.closer.hex_day1_ratioCorrect,0.5,'Tail','right')
[~, pCloserHexDay2] = ttest(D.closer.hex_day2_ratioCorrect,0.5,'Tail','right')
[~, pCloserClusDay3] = ttest(D.closer.clus_day3_ratioCorrect,0.5,'Tail','right')
[~, pCloserClusDay4] = ttest(D.closer.clus_day4_ratioCorrect,0.5,'Tail','right')

% stats for improvement beween days
[~, pCloserHexDays] = ttest2(D.closer.hex_day2_ratioCorrect,D.closer.hex_day1_ratioCorrect,'Tail','right')
[~, pCloserClusDays] = ttest2(D.closer.clus_day4_ratioCorrect,D.closer.clus_day3_ratioCorrect,'Tail','right')


%% navigation task


for iDist = 1:3
    subplot(2,3,iDist+3)
    hold on
    
    title(['Navig dist ' num2str(iDist + 1)]) % between 2 and 4
    
    % individual datapoints for current initial distance
    day1 = mean(D.navig.hex_day1_nSteps(:,:,iDist),2);
    day2 = mean(D.navig.hex_day2_nSteps(:,:,iDist),2);
    day3 = mean(D.navig.clus_day3_nSteps(:,:,iDist),2);
    day4 = mean(D.navig.clus_day4_nSteps(:,:,iDist),2);
    
    % Plot individual data points
    scatter(ones(nSub,1),day1,3,'filled','r','MarkerFaceAlpha',0.5)
    scatter(2*ones(nSub,1),day2,3,'filled','r','MarkerFaceAlpha',0.5)
    scatter(3*ones(nSub,1),day3,3,'filled','b','MarkerFaceAlpha',0.5)
    scatter(4*ones(nSub,1),day4,3,'filled','b','MarkerFaceAlpha',0.5)
    
    % plot means and SEM error bars
    wBar = 0.15;
    plot([1+0.2-wBar, 1+0.2+wBar], [mean(day1),mean(day1)],'k','lineWidth',2)
    plot([2-0.2-wBar, 2-0.2+wBar], [mean(day2),mean(day2)],'k','lineWidth',2)
    plot([3+0.2-wBar, 3+0.2+wBar], [mean(day3),mean(day3)],'k','lineWidth',2)
    plot([4-0.2-wBar, 4-0.2+wBar], [mean(day4),mean(day4)],'k','lineWidth',2)
    
    errorbar(1+0.2,mean(day1),std(day1)/sqrt(nSub),'k')
    errorbar(2-0.2,mean(day2),std(day2)/sqrt(nSub),'k')
    errorbar(3+0.2,mean(day3),std(day3)/sqrt(nSub),'k')
    errorbar(4-0.2,mean(day4),std(day4)/sqrt(nSub),'k')
    
    % Plot lines connecting the data points
    for iSub=1:nSub
        plot([1,2],[day1(iSub),day2(iSub)],'color',[1 0 0 0.2])
        plot([3,4],[day3(iSub),day4(iSub)],'color',[0 0 1 0.2])
    end
    
    xlim([0,5]);
    ylabel('nSteps')
    labels = {'hex day1','hex day2','comm day3', 'comm day4'};
    set(gca,'FontSize',14,'XTick',1:4,'XTickLabels',labels,'XTickLabelRotation', 45)
    
    % stats for whether answers oneach day were different from chanve
    [~, pNavigDay1] = ttest(day1,0.5,'Tail','right');
    sprintf('p initDist %d, day 1 = %d',iDist+1,pNavigDay1)
    [~, pNavigDay2] = ttest(day2,0.5,'Tail','right');
    sprintf('p initDist %d, day 2 = %d',iDist+1,pNavigDay2)
    [~, pNavigDay3] = ttest(day3,0.5,'Tail','right');
    sprintf('p initDist %d, day 3 = %d',iDist+1,pNavigDay3)
    [~, pNavigDay4] = ttest(day4,0.5,'Tail','right');
    sprintf('p initDist %d, day 4 = %d',iDist+1,pNavigDay4)
    
    % stats for improvement beween days
    [~, pNavigHexDays] = ttest2(day1,day2,'Tail','right');
    [~, pNavigClusDays] = ttest2(day3,day4,'Tail','right');
    sprintf('p initDist %d, diff days Hex = %d',iDist+1,pNavigHexDays)
    sprintf('p initDist %d, diff days Comm = %d',iDist+1,pNavigClusDays)
end






