function [mean_coh] = normCohByPairings(all_coh, targetIdxs,referenceIdxs)
%NORMCOHBYPAIRINGS Summary of this function goes here
%   Detailed explanation goes here

[targetPairedIdxs,referencePairedIdxs] = ndgrid(targetIdxs,referenceIdxs);
Z = [targetPairedIdxs(:),referencePairedIdxs(:)];
out_coh = nan(size(all_coh,1), size(all_coh,2), size(Z,1));
for zz = 1: size(Z,1)
    
    target_coh = all_coh(:, :, targetPairedIdxs(zz), :);
    reference_coh = all_coh(:, :, referencePairedIdxs(zz), :);
    
    norm_coh = target_coh-reference_coh;
    
    out_coh(:, :, zz) = norm_coh;
    
end

mean_coh = nanmean(out_coh, 3);
end

