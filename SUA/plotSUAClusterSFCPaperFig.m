clearvars
clc

%% Paths to COH matrices
data_root = '/Volumes/external3/data/changlab/ptp-vicarious-reward/sfcoh-split-by-cell-type';
path2CohMat  = fullfile( data_root, 'all-trialtype_K2_SFC_Matrix.mat' );
figStrAdd = 'Cluster2';
specialStr = 'trialtype';


%% Settings

cnorm = true;

load(path2CohMat);


%% Post-load settings
minPlotFreq = 10;
maxPlotFreq = 60;
minPlotMs = -150;
maxPlotMs = 600;
smoothFac = 3;
ROI_MinFreq = 30; % 30-50 or 8-20
ROI_MaxFreq = 50;
ROI_MinFreq = 10; % 30-50 or 8-20
ROI_MaxFreq = 20;
ROI_MinMs = 50;
ROI_MaxMs = 300;
saveFigs = true;
useMedian4Scatter = false;

%% Reduce size of averageMat
validFreqIdxs = find(loadStruct.coh_f >= minPlotFreq & loadStruct.coh_f < maxPlotFreq);
meanTimeBins = cellfun(@mean,loadStruct.binned_t)+loadStruct.t(1);
validTimeIdxs = find(meanTimeBins > minPlotMs & meanTimeBins < maxPlotMs);
SFCMAT = SFCMAT(validFreqIdxs, validTimeIdxs, :, :, :);
timeLabels = meanTimeBins(validTimeIdxs);
freqLabels = loadStruct.coh_f(validFreqIdxs);

%% Split apart into conditions and calculate mean
t1_cued_normSelf_popCoh         = squeeze(SFCMAT(:,:,1,1,:));
t1_cued_normOther_popCoh        = squeeze(SFCMAT(:,:,1,2,:));
t1_choice_normSelf_popCoh       = squeeze(SFCMAT(:,:,1,3,:));
t1_choice_normOther_popCoh      = squeeze(SFCMAT(:,:,1,4,:));

t2_cued_normSelf_popCoh         = squeeze(SFCMAT(:,:,2,1,:));
t2_cued_normOther_popCoh        = squeeze(SFCMAT(:,:,2,2,:));
t2_choice_normSelf_popCoh       = squeeze(SFCMAT(:,:,2,3,:));
t2_choice_normOther_popCoh      = squeeze(SFCMAT(:,:,2,4,:));

t1_cued_normSelf_mCoh         = nanmean(t1_cued_normSelf_popCoh,3);
t1_cued_normOther_mCoh        = nanmean(t1_cued_normOther_popCoh,3);
t1_choice_normSelf_mCoh       = nanmean(t1_choice_normSelf_popCoh,3);
t1_choice_normOther_mCoh      = nanmean(t1_choice_normOther_popCoh,3);

t2_cued_normSelf_mCoh         = nanmean(t2_cued_normSelf_popCoh,3);
t2_cued_normOther_mCoh        = nanmean(t2_cued_normOther_popCoh,3);
t2_choice_normSelf_mCoh       = nanmean(t2_choice_normSelf_popCoh,3);
t2_choice_normOther_mCoh      = nanmean(t2_choice_normOther_popCoh,3);

%% Plot Figure 2 ACC --> BLA

fig2H = figure(2);
clf

ROI_FreqIdxs = find(freqLabels > ROI_MinFreq & freqLabels < ROI_MaxFreq);
ROI_TimeIdxs= find(timeLabels > ROI_MinMs & timeLabels < ROI_MaxMs);
ROI_FreqRange = freqLabels(ROI_FreqIdxs(end))-freqLabels(ROI_FreqIdxs(1));
ROI_TimeRange = timeLabels(ROI_TimeIdxs(end))-timeLabels(ROI_TimeIdxs(1));



axH(3) = subplot(3,3,3);
hold on;
t1_choice_normSelf_ROICoh = t1_choice_normSelf_popCoh(ROI_FreqIdxs,:,:);
t1_choice_normSelf_ROICoh_avgXFrq = squeeze(nanmean(t1_choice_normSelf_popCoh(ROI_FreqIdxs,:,:),1));
t1_choice_normOther_ROICoh = t1_choice_normOther_popCoh(ROI_FreqIdxs,:,:);
t1_choice_normOther_ROICoh_avgXFrq = squeeze(nanmean(t1_choice_normOther_popCoh(ROI_FreqIdxs,:,:),1));
[choice_self_lineOut, ~] = stdshade(t1_choice_normSelf_ROICoh_avgXFrq',.5,'b', timeLabels, smoothFac);
[choice_other_lineOut, ~] = stdshade(t1_choice_normOther_ROICoh_avgXFrq',.5,'r', timeLabels, smoothFac);
axis tight
hline(0, 'k:');
vline(0, 'k-');
ylim([-.01 0.01])
legend([choice_self_lineOut choice_other_lineOut],{'Self','Other'})
title(sprintf('%3.1f to %3.1f Hz', ROI_MinFreq, ROI_MaxFreq));





axH(6) = subplot(3,3,6);
hold on
t1_cued_normSelf_ROICoh = t1_cued_normSelf_popCoh(ROI_FreqIdxs,:,:);
t1_cued_normSelf_ROICoh_avgXFrq = squeeze(nanmean(t1_cued_normSelf_popCoh(ROI_FreqIdxs,:,:),1));
t1_cued_normOther_ROICoh = t1_cued_normOther_popCoh(ROI_FreqIdxs,:,:);
t1_cued_normOther_ROICoh_avgXFrq = squeeze(nanmean(t1_cued_normOther_popCoh(ROI_FreqIdxs,:,:),1));
[cued_self_lineOut, ~] = stdshade(t1_cued_normSelf_ROICoh_avgXFrq',.5,'b', timeLabels, smoothFac);
[cued_other_lineOut, ~] = stdshade(t1_cued_normOther_ROICoh_avgXFrq',.5,'r', timeLabels, smoothFac);
axis tight
hline(0, 'k:');
vline(0, 'k-');
ylim([-.01 0.01])
legend([cued_self_lineOut cued_other_lineOut],{'Self (Forced) - Bottle(Forced)','Other-Bottle'})
title(sprintf('%3.1f to %3.1f Hz', ROI_MinFreq, ROI_MaxFreq));



axH(7) = subplot(3,3,7);
hold on;

choice_self_ROI_vals  = squeeze(nanmean(nanmean(t1_choice_normSelf_popCoh(ROI_FreqIdxs, ROI_TimeIdxs, :),1),2));
choice_self_ROI_vals(isnan(choice_self_ROI_vals)) = [];
choice_self_ROI_mean = nanmean(choice_self_ROI_vals);
choice_self_ROI_var = nanstd(choice_self_ROI_vals)/sqrt(numel(~isnan(choice_self_ROI_vals)));
bar(1, choice_self_ROI_mean, 'FaceColor', 'b', 'EdgeColor', 'none');
errorbar(1,choice_self_ROI_mean,choice_self_ROI_var, 'k', 'LineWidth', 2)

choice_other_ROI_vals  = squeeze(nanmean(nanmean(t1_choice_normOther_popCoh(ROI_FreqIdxs, ROI_TimeIdxs, :),1),2));
choice_other_ROI_vals(isnan(choice_other_ROI_vals)) = [];
choice_other_ROI_mean = nanmean(choice_other_ROI_vals);
choice_other_ROI_var = nanstd(choice_other_ROI_vals)/sqrt(numel(~isnan(choice_other_ROI_vals)));
bar(2, choice_other_ROI_mean, 'FaceColor', 'r', 'EdgeColor', 'none');
errorbar(2, choice_other_ROI_mean,choice_other_ROI_var, 'k', 'LineWidth', 2)


xticks([1 2]);
xticklabels({'Self-Bottle', 'Other-Bottle'})
xtickangle(45)
[p,~,~] = ranksum(choice_self_ROI_vals,choice_other_ROI_vals);
title(sprintf('Wilcoxon rank sum p=%6.5f', p));



axH(8) = subplot(3,3,8);
hold on;


cued_self_ROI_vals = squeeze(nanmean(nanmean(t1_cued_normSelf_popCoh(ROI_FreqIdxs, ROI_TimeIdxs, :),1),2));
cued_self_ROI_vals(isnan(cued_self_ROI_vals)) = [];
cued_self_ROI_mean = nanmean(cued_self_ROI_vals);
cued_self_ROI_var = nanstd(cued_self_ROI_vals)/sqrt(numel(~isnan(cued_self_ROI_vals)));
bar(1, cued_self_ROI_mean, 'FaceColor', 'b', 'EdgeColor', 'none');
errorbar(1,cued_self_ROI_mean,cued_self_ROI_var, 'k', 'LineWidth', 2)

cued_other_ROI_vals = squeeze(nanmean(nanmean(t1_cued_normOther_popCoh(ROI_FreqIdxs, ROI_TimeIdxs, :),1),2));
cued_other_ROI_vals(isnan(cued_other_ROI_vals)) = [];
cued_other_ROI_mean = nanmean(cued_other_ROI_vals);
cued_other_ROI_var = nanstd(cued_other_ROI_vals)/sqrt(numel(~isnan(cued_other_ROI_vals)));
bar(2, cued_other_ROI_mean, 'FaceColor', 'r', 'EdgeColor', 'none');
errorbar(2,cued_other_ROI_mean,cued_other_ROI_var, 'k', 'LineWidth', 2)

xticks([1 2]);
xticklabels({'Self (Forced)', 'Other (Forced)'})
xtickangle(45)
[p,~,~] = ranksum(cued_self_ROI_vals,cued_other_ROI_vals);
title(sprintf('Wilcoxon rank sum p=%6.5f', p));



titleStr = sprintf('ACC Spikes-BLA FP - %3.1f to %3.1f Hz - %4d ms to %4d ms', ROI_MinFreq, ROI_MaxFreq, ROI_MinMs, ROI_MaxMs);
sgtitle(titleStr, 'FontSize', 14)
% 
% fig2saveStr = strrep(sprintf('Fig2_ACCSpikes-BLAFP_%3.1f_%3.1fHz_%d-%d_%s_Z_%s', ROI_MinFreq, ROI_MaxFreq,ROI_MinMs, ROI_MaxMs, specialStr, figStrAdd),'.','_');
% fig2saveStr = [fig2saveStr, '.pdf'];
% if saveFigs
%     fig2SavePath = fullfile(pwd, fig2saveStr);
%     makePortraitPDF(fig2H, fig2SavePath);
% end





