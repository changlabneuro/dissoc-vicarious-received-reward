function [cMin, cMax] = normSubsetCAxes(subsetOfAxes)
%NORMSUBSETCAXES Summary of this function goes here

CLims = nan(size(subsetOfAxes,1), 2);

for ax = 1:size(subsetOfAxes,1)
    
    thisAx = subsetOfAxes(ax);
   
    CLims(ax,:) = thisAx.CLim;
end

cMin = min(CLims(:,1));
cMax = max(CLims(:,2));


for ax = 1:size(subsetOfAxes,1)
    
    thisAx = subsetOfAxes(ax);
   
    caxis(thisAx,[cMin cMax])
end
