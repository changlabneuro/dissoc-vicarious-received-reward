function [aHand,bHand,cHand, aVals, bVals,cVals] = plotTripartSpectroTimeSeriesROI(axeH, aMat,bMat, cMat, freqIdxs, timeLabels)
%PLOTSPECTROTIMESERIESROI Summary of this function goes here
%   Detailed explanation goes here


smoothFac = 3;
axes(axeH);


%%
aMat_avgXFrq =squeeze(nanmean(real(aMat(:,freqIdxs,:)),2));
bMat_avgXFrq = squeeze(nanmean(real(bMat(:,freqIdxs,:)),2));
cMat_avgXFrq = squeeze(nanmean(real(cMat(:,freqIdxs,:)),2));

hold on
[aHand, ~] = stdshade(aMat_avgXFrq,.5,'r', timeLabels, smoothFac);
[bHand, ~] = stdshade(bMat_avgXFrq,.5,'b', timeLabels, smoothFac);
[cHand, ~] = stdshade(cMat_avgXFrq,.5,'k', timeLabels, smoothFac);
axis tight
hline(0, 'k:');
vline(0, 'k-');

%%

aVals = aMat_avgXFrq;
bVals = bMat_avgXFrq;
cVals = cMat_avgXFrq;
end

