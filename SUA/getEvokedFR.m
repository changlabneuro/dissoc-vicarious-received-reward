function [binnedHz, celledHz] = getEvokedFR(spikeS,eventS, binS, windowS)
%GETEVOKEDFR Summary of this function goes here
%   Detailed explanation goes here
% binS(1,:) = binStartS;
% binS(2,:) = binStopS;
binStartS = binS(1,:);
binStopS = binS(2,:);
numBins = size(binS,2);
binSizeS = mode(binStopS - binStartS); % assume uniform bins....
celledHz = cell(numel(eventS), 0);
        event_binC = nan(numel(eventS), numBins);
        for e = 1:numel(eventS)
            thisEventS = eventS(e);
            windowStartS = thisEventS - (windowS/2);
            windowStopS = thisEventS + (windowS/2);
            windowSpikeIdxs = spikeS >windowStartS & spikeS < windowStopS;
            windowSpikeAbsS = spikeS(windowSpikeIdxs);
            windowSpikeRelS = (windowSpikeAbsS - windowStartS) - (windowS/2) ;
            celledHz{e,1} = windowSpikeRelS';
            % --------------- Spit up this trial spikes into bins ---------------
            for b = 1:numBins
                thisBinStartS = binStartS(b);
                thisBinsStopS = binStopS(b);
                thisBinSpikeS = windowSpikeRelS(windowSpikeRelS >= thisBinStartS & windowSpikeRelS < thisBinsStopS);
                event_binC(e, b) = numel(thisBinSpikeS);
            end
        end
        binnedHz = event_binC/binSizeS;
end