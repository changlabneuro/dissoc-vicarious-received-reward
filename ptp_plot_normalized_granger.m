%%
clear all
clc
 
load('GrangerPrePro.mat');
 
 
ROI_MinFreq = 30;
ROI_MaxFreq = 50;
%%
 
 
norm_self_acc_bla_mat = self_acc_bla_mat-none_acc_bla_mat;
norm_self_bla_acc_mat = self_bla_acc_mat-none_bla_acc_mat;
 
norm_other_acc_bla_mat = other_acc_bla_mat-none_acc_bla_mat;
norm_other_bla_acc_mat = other_bla_acc_mat-none_bla_acc_mat;
 
tsMin = 0;
tsMax = 0.1;
 
ROI_TimeIdxs = find(tLabs >= ROI_MinMs & tLabs <= ROI_MaxMs);
ROI_FreqIdxs = find(fLabs >= ROI_MinFreq & fLabs <= ROI_MaxFreq);
 
ROI_TimeRange = tLabs(ROI_TimeIdxs(end))-tLabs(ROI_TimeIdxs(1));
ROI_FreqRange = fLabs(ROI_FreqIdxs(end))-fLabs(ROI_FreqIdxs(1));
 
figH = figure(1);
clf
 
axesH(1) = subplot(2,1,1);
title({'Self-None vs Other-None trials', 'ACC ---> BLA'});
[self_acc_bla_H,other_acc_bla_H, self_acc_bla_Vals, other_acc_bla_Vals] = plotSpectroTimeSeriesROI(axesH(1) , norm_self_acc_bla_mat,norm_other_acc_bla_mat, ROI_FreqIdxs, tLabs);
xlim([-500 500]);
xlabel('Time (ms)');
legend([self_acc_bla_H other_acc_bla_H],{'Self-None','Other-none'},'Location', 'best')
ylim([tsMin tsMax]);
ylim([-0.05 .09])
vline(0, 'k-')
axesH(2) = subplot(2,1,2);
 
title({'Self-None vs Other-None trials', 'BLA ---> ACC'});
[self_bla_acc_H,other_bla_acc_H, self_bla_acc_Vals, other_bla_acc_Vals] = plotSpectroTimeSeriesROI(axesH(2) , norm_self_bla_acc_mat,norm_other_bla_acc_mat, ROI_FreqIdxs, tLabs);
xlim([-500 500]);
xlabel('Time (ms)');
legend([self_bla_acc_H other_bla_acc_H],{'Self-None','Other-none'},'Location', 'best')
ylim([tsMin tsMax]);
ylim([-0.05 .09])
vline(0, 'k-')
makeLandscapePDF(figH, 'LowGammaGrangerDirectionalityNormalizedByBottle.pdf')