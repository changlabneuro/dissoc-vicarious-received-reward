function [evokedFRs, binMeanS] = getEvokedFRs(spikeTs,eventTs, conditionIdxs, windowMs, binSizesMs)
%GETCONDITIONALFRS Summary of this function goes here
%  Detailed explanation goes here
 
numConditions = size(conditionIdxs,2);
paddingFac = 0.3;
windowS = windowMs/1e3;
binsEarlierTimeS      = -(windowS/2) + (-(windowS/2) *paddingFac);
binsLaterTimeS       = (windowS/2) + ((windowS/2)*paddingFac);
binSizeMs          = binSizesMs(1);
binSlideMs         = binSizesMs(2);
binSlideS          = binSlideMs/1e3;
binSizeS          = binSizeMs/1e3;
binStartS          = (binsEarlierTimeS) : binSlideS : binsLaterTimeS;
numBins           = size(binStartS,2);
binStopS          = nan(size(binStartS));
binMeanS          = nan(size(binStartS));
 
for b = 1 : numBins
  binStopS(b) = binStartS(b) + binSizeS;
  binMeanS(b) = mean([binStartS(b) binStopS(b)]);
end
 
evokedFRs = cell(numConditions,1);
 
for c = 1:numConditions
   
  cIdxs = conditionIdxs{c};
   
  cEventTs = eventTs(cIdxs);
   
  eventSpikeS = cell(numel(cEventTs),1);
   
  eventBinnedCounts = nan(numel(cEventTs), numBins);
   
  for e = 1:numel(cEventTs)
     
    eventTimeS = cEventTs(e);
     
    eventWindowStartS = (eventTimeS - ( (windowS/2) + ((windowS/2)*paddingFac)) );
     
    eventWindowStopS = (eventTimeS + ( (windowS/2) + ((windowS/2)*paddingFac)) ) + binSizeS;
     
    windowSpikeRelS = spikeTs(spikeTs >= eventWindowStartS & spikeTs <= eventWindowStopS) - (eventTimeS);
     
    eventSpikeS{e} = windowSpikeRelS;
     
    % --------------- Spit up this trial spikes into bins ---------------
    for b = 1:numBins
       
      thisBinStartS = binStartS(b);
      thisBinsStopS = binStopS(b);
      thisBinSpikeS = windowSpikeRelS(windowSpikeRelS >= thisBinStartS & windowSpikeRelS < thisBinsStopS);
      eventBinnedCounts(e, b) = numel(thisBinSpikeS);
    end
     
  end
   
  eventBinnedHz = eventBinnedCounts/binSizeS;
   
  evokedFRs{c}(1:numBins, 1:size(eventBinnedHz,1)) = eventBinnedHz';
   
   
 
end
 
end