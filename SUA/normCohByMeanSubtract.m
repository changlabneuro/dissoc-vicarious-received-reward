function [mean_coh] = normCohByMeanSubtract(all_coh, targetIdxs,referenceIdxs)
%NORMCOHBYPAIRINGS Summary of this function goes here
%   Detailed explanation goes here

target_coh = all_coh(:,:,targetIdxs);

ref_coh = all_coh(:,:,referenceIdxs);

mean_target_coh = nanmean(target_coh,3);
mean_ref_coh = nanmean(ref_coh,3);

mean_coh =mean_target_coh-mean_ref_coh;
end

