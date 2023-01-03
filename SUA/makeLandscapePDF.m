function makeLandscapePDF(figureHandle,savePath)
%MAKELANDSCAPEPDF Summary of this function goes here
%   Detailed explanation goes here

set(figureHandle,'PaperOrientation','landscape');
set(figureHandle,'PaperUnits','normalized');
set(figureHandle,'PaperPosition', [0 0 1 1]);
print(figureHandle, '-dpdf', savePath);

end

