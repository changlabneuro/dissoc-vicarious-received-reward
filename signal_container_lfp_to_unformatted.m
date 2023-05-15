function [data, labels, categories, t] = signal_container_lfp_to_unformatted(lfp_container)

%   SIGNAL_CONTAINER_LFP_TO_UNFORMATTED 
%
%     [data, labels, categories, t] = ...
%       signal_container_lfp_to_unformatted( lfp_container );
%
%     decomposes the SignalContainer object `lfp_container` into its
%     constituent parts, using only builtin Matlab data types. 
%
%     . `data` is an MxN double array of M trials by N samples @ 1khz.
%     . `labels` is an MxN' categorical array of M trials by N' categories,
%         identifying rows of `data`.
%     . `categories` is a N'x1 cell array of strings identifying the 
%         columns of `labels`.
%     . `t` is a 1xN double vector giving the time in ms relative to the 
%       source epoch onset, for each sample of `data`.
%
%     This function depends on the `dsp`, `global`, and `categorical`
%     repositories, available on github:
%
%       . dsp:          https://github.com/nfagan/dsp
%       . global:       https://github.com/nfagan/global
%       . categorical:  https://github.com/nfagan/categorical
%     
%     See also SignalContainer, SparseLabels, fcat

[labels, categories] = categorical( fcat.from(lfp_container.labels) );
data = lfp_container.data;

t_props = lfp_container.get_time_props();
t_min = t_props(1);
t = t_min:(t_min + size(data, 2)-1);

end