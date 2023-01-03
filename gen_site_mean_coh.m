file_list = shared_utils.io.findmat( '/Volumes/external3/data/changlab/ptp-vicarious-reward/newCoh2' );

% match = contains( file_list, 'acc_unit_index') & contains( file_list, 'bla_COH' );
% file_list = file_list(match);

%%

mean_each = { 'outcome', 'trial-type', 'administration' };
mask_func = @(l) findnone( l, 'errors' );

do_save = true;
dst_p = '/Volumes/external3/data/changlab/ptp-vicarious-reward/site-mean-sfcoh';

for i = 1:numel(file_list)
  
shared_utils.general.progress( i, numel(file_list) );
  
dst_filename = shared_utils.io.filenames( file_list{i}, true );

coh = load( file_list{i} );
coh_mat = permute( coh.cohMatrix, [3, 1, 2] );
trial_labs = make_trial_labels( coh.ccLFPInfo, coh.ccSpikeInfo );

f = coh.coh_f;
f_ind = f < 100;
f = f(f_ind);
coh_mat = coh_mat(:, f_ind, :);

bin_t = shared_utils.vector.slidebin( coh.t, 150, 50, true );
assert( numel(bin_t) == size(coh_mat, 3) );
t = cellfun( @median, bin_t );

mask = mask_func( trial_labs );
[mean_labs, I] = keepeach( trial_labs', mean_each, mask );
means = rowop( coh_mat, I, @(x) nanmean(x, 1) );

assert_ispair( means, mean_labs );
assert( size(means, 2) == numel(f) && size(means, 3) == numel(t) );

to_save = struct( ...
    'coh', means ...
  , 'labels', mean_labs ...
  , 'f', f ...
  , 't', t ...
);

save_p = fullfile( dst_p, dst_filename );
if ( do_save )
  fprintf( '\n Saving "%s".', save_p );
  save( save_p, 'to_save' );
end

end

%%

function trial_labs = make_trial_labels(lfp_info, spk_info)

trial_labs = fcat.from( lfp_info ...
  , {'day', 'lfp-channel', 'magnitude', 'outcome' ...
  , 'trial-type', 'administration', 'lfp-region'} );
replace( trial_labs, 'acc', 'lfp-acc' );
replace( trial_labs, 'bla', 'lfp-bla' );

sf_labs = fcat.from( spk_info ...
  , {'day', 'spk-channel', 'spk-region', 'unit-index'} );
replace( sf_labs, 'acc', 'spk-acc' );
replace( sf_labs, 'bla', 'spk-bla' );

join( trial_labs, sf_labs );

end
