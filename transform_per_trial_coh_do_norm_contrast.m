function [coh, labels, f, t] = transform_per_trial_coh_do_norm_contrast(coh, labels, f, t, each, across, varargin)

[coh, labels, f, t] = transform_per_trial_coh( coh, labels, f, t );
[coh, labels] = do_norm_contrast( coh, labels, findall(labels, each), across, varargin{:} );

end
