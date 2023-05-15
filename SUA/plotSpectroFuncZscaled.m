function plotSpectroFuncZscaled(orig_Mat, axesH, timeLabels,freqLabels)
%PLOTSPECTROFUNC Summary of this function goes here
%   Detailed explanation goes here
axes(axesH)
colormap jet

reverseOrder = false;

if reverseOrder
    f_Mat = imgaussfilt(orig_Mat,2);
    z_Mat = zscore(f_Mat);
    p_Mat = rescale(z_Mat,min(min(orig_Mat)),max(max(orig_Mat)));
else
    z_Mat = zscore(orig_Mat);
    f_Mat = imgaussfilt(z_Mat,2);
    p_Mat = rescale(f_Mat,min(min(orig_Mat)),max(max(orig_Mat)));
end

avoid_zscore = true;  % 03/09/22
if ( avoid_zscore )
  p_Mat = imgaussfilt( orig_Mat, 2 );
end


imagesc('XData',timeLabels,'YData',freqLabels,'CData',p_Mat)
set(gca,'YDir','normal')
axis tight
vline(0, 'k-')
colorbar


end

