function [aHand,bHand, aVals, bVals] = plotSpectroTimeSeriesROI(axeH, aMat,bMat, freqIdxs, timeLabels)
%PLOTSPECTROTIMESERIESROI Summary of this function goes here
%   Detailed explanation goes here


smoothFac = 3;
axes(axeH);


%%
aMat_avgXFrq =squeeze(nanmean(real(aMat(:,freqIdxs,:)),2));
bMat_avgXFrq = squeeze(nanmean(real(bMat(:,freqIdxs,:)),2));

hold on
[aHand, ~] = stdshade(aMat_avgXFrq,.5,'b', timeLabels, smoothFac);
[bHand, ~] = stdshade(bMat_avgXFrq,.5,'r', timeLabels, smoothFac);
axis tight
hline(0, 'k:');
vline(0, 'k-');

%%

aVals = aMat_avgXFrq;
bVals = bMat_avgXFrq;
end

