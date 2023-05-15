sub_ps = strrep( subPaths, '\', '/' );
match_fs = shared_utils.io.filenames( sub_ps, true );

search_p = '/Volumes/external3/data/changlab/ptp-vicarious-reward/newCoh2';
exist_files = false( size(match_fs) );
sub_fs = cell( size(match_fs) );

for i = 1:numel(match_fs)
  search_f = fullfile( search_p, match_fs{i} );
  exist_files(i) = shared_utils.io.fexists( search_f );
  sub_fs{i} = search_f;
end


%%

normMethod = 'STDNORM';
sameRegionCount = 0;
vc = 0;
t1c = 0;
t2c = 0;

site_info = {};

for i = 1:size(sub_fs,1)
  
fprintf( '\n%d of %d', i, size(sub_fs, 1) );

loadStruct = load(sub_fs{i});

if ~strcmp(loadStruct.spkReg, loadStruct.lfpRegs)

    vc = vc+1;

    if strcmp(loadStruct.spkReg, 'acc') && strcmp(loadStruct.lfpRegs, 'bla')
        t1c = t1c+1;
        type = 1;
    elseif strcmp(loadStruct.spkReg, 'bla') && strcmp(loadStruct.lfpRegs, 'acc')
        t2c = t2c+1;
        type = 2;
    else
        error('bad, not rad')
        type = [];
    end
    
    %%  just get pair info
    copy_fields = { 'lfpChans', 'lfpDay', 'lfpRegs', 'spikeIdx', 'spkChan', 'spkReg' };
    field_subset = cellfun( @(x) loadStruct.(x), copy_fields, 'un', 0 );
    site_info(end+1, :) = field_subset;
%     if ( 1 ), continue; end
    
    %%

    all_coh  = loadStruct.cohMatrix;
    administration = loadStruct.ccLFPInfo(:,6);
    outcomes = loadStruct.ccLFPInfo(:,4);
    trialtypes = loadStruct.ccLFPInfo(:,5);
    
%     trial_ind = find( trial_labels, loadStruct.lfpDay );
%     assert( numel(trial_ind) == size(loadStruct.cohMatrix, 3) );
    
    
    preIdxs = find(strcmp(administration, 'pre'));
    bothIdxs = find(strcmp(outcomes, 'both'));
    selfIdxs = find(strcmp(outcomes, 'self'));
    otherIdxs = find(strcmp(outcomes, 'other'));
    noneIdxs = find(strcmp(outcomes, 'none'));
    cuedIdxs = intersect(find(strcmp(trialtypes, 'cued')), preIdxs);
    choiceIdxs =  intersect(find(strcmp(trialtypes, 'choice')), preIdxs);

    % Actual indexs we're interested in
    cued_both = intersect(cuedIdxs, bothIdxs);
    cued_self = intersect(cuedIdxs, selfIdxs);
    cued_other = intersect(cuedIdxs, otherIdxs);
    cued_none = intersect(cuedIdxs, noneIdxs);
    choice_both = intersect(choiceIdxs, bothIdxs);
    choice_self = intersect(choiceIdxs, selfIdxs);
    choice_other = intersect(choiceIdxs, otherIdxs);
    choice_none = intersect(choiceIdxs, noneIdxs);
    
    if ( 1 )
      i_sn = [ cued_self(:); cued_none(:) ];
      mu_s = nanmean( all_coh(:, :, i_sn), 3 );
      mu_nn = mu_s;

      i_on = [ cued_other(:); cued_none(:) ];
      mu_o = nanmean( all_coh(:, :, i_on), 3 );
      mu_n = mu_o;
    else
      mu_s = nanmean( all_coh(:, :, cued_self), 3 );
      mu_n = nanmean( all_coh(:, :, cued_none), 3 );
      mu_nn = mu_n;
      mu_o = nanmean( all_coh(:, :, cued_other), 3 );
    end
    
    targ_subset_sn = src.averageMat(:, :, 1, 1, vc);
    targ_subset_on = src.averageMat(:, :, 1, 2, vc);
    delta_sn = (nanmean(all_coh(:, :, cued_self), 3) - mu_s) - ...
      (nanmean(all_coh(:, :, cued_none), 3) - mu_nn);
    delta_on = (nanmean(all_coh(:, :, cued_other), 3) - mu_o) - ...
      (nanmean(all_coh(:, :, cued_none), 3) - mu_n);

    d_sn = max( columnize(abs(delta_sn - targ_subset_sn)) );
    d_on = max( columnize(abs(delta_on - targ_subset_on)) );
    assert( ...
      (all(isnan(targ_subset_sn(:))) && all(isnan(targ_subset_on(:)))) || ...
      (all(targ_subset_sn(:) == 0) && all(targ_subset_on(:) == 0)) || ...
      (d_sn < 1e-4 && d_on < 1e-4) );
    
%     if strcmp(normMethod, 'STDNORM')
%         averageMat(:, :, type, 1, vc)  = normCohByMeanSubtract(all_coh, cued_self,  cued_none);
%         averageMat(:, :, type, 2, vc)  = normCohByMeanSubtract(all_coh, cued_other, cued_none);
%         averageMat(:, :, type, 3, vc)  = normCohByMeanSubtract(all_coh, choice_self,  choice_none);
%         averageMat(:, :, type, 4, vc)  = normCohByMeanSubtract(all_coh, choice_other, choice_none);
% 
%     elseif strcmp(normMethod, 'PAIRNORM')
%         averageMat(:, :, type, 1, vc)  = normCohByPairings(all_coh, cued_self,  cued_none);
%         averageMat(:, :, type, 2, vc)  = normCohByPairings(all_coh, cued_other, cued_none);
%         averageMat(:, :, type, 3, vc)  = normCohByPairings(all_coh, choice_self,  choice_none);
%         averageMat(:, :, type, 4, vc)  = normCohByPairings(all_coh, choice_other, choice_none);
%     else
%         error('Select a valid normalization method');
%     end
else
    sameRegionCount = sameRegionCount+1;
    sameRegionPaths{sameRegionCount, 1} = sub_fs{i};
end

end

%%  split by social gaze

cons = dsp3.get_consolidated_data( dsp3.set_dataroot('/Volumes/external3/data/changlab/dsp3') );
event_labels = fcat.from( cons.events.labels );
trial_labels = fcat.from( cons.trial_data.labels );
% @NOTE: Bottle and monkey are switched in raw data.
late_monk_look_counts = cons.trial_data.data(:, cons.trial_key('lateBottleLookCount'));
late_bottle_look_counts = cons.trial_data.data(:, cons.trial_key('lateLookCount'));

monk_inds = find( late_monk_look_counts > 0 );
bottle_inds = find( late_bottle_look_counts > 0 );
both_inds = intersect( monk_inds, bottle_inds );
monk_inds = setdiff( monk_inds, both_inds );
bottle_inds = setdiff( bottle_inds, both_inds );

num_saved = 0;

for i = 1:size(sub_fs, 1)
  
fprintf( '\n%d of %d', i, size(sub_fs, 1) );

loadStruct = load(sub_fs{i});

if strcmp(loadStruct.spkReg, loadStruct.lfpRegs)
  continue;
end

trial_ind = find( trial_labels, loadStruct.lfpDay );
assert( numel(trial_ind) == size(loadStruct.cohMatrix, 3) );

[~, i_bottle] = intersect( trial_ind, bottle_inds, 'stable' );
[~, i_monk] = intersect( trial_ind, monk_inds, 'stable' );
[~, i_both] = intersect( trial_ind, both_inds, 'stable' );

coh_labs = make_sf_labels( loadStruct );
addcat( coh_labs, 'looks_to' );
setcat( coh_labs, 'looks_to', 'looks-to-bottle', i_bottle );
setcat( coh_labs, 'looks_to', 'looks-to-monkey', i_monk );
setcat( coh_labs, 'looks_to', 'looks-to-bottle-and-monkey', i_both );

coh_mat = permute( loadStruct.cohMatrix, [3, 1, 2] );
sub_filename = shared_utils.io.filenames( sub_fs{i}, true );
sub_p = fullfile( fileparts(fileparts(sub_fs{i})), 'per-trial-sfcoh', sub_filename );

keep_f = loadStruct.coh_f < 100;
coh_mat = coh_mat(:, keep_f, :);
coh_f = loadStruct.coh_f(keep_f);

coh_t = shared_utils.vector.slidebin( loadStruct.t, 150, 50, true );
coh_t = cellfun( @median, coh_t );
assert( numel(coh_t) == size(coh_mat, 3) );
assert( numel(coh_f) == size(coh_mat, 2) );

num_saved = num_saved + 1;
fprintf( ' Saving %d', num_saved );
save( sub_p, 'coh_mat', 'coh_f', 'coh_labs', 'coh_t' );

end

%%

norm_across = { 'magnitudes' };

load_mats = shared_utils.io.findmat( '/Volumes/external3/data/changlab/ptp-vicarious-reward/per-trial-sfcoh' );
[coh, coh_labels, f, t] = bfw.load_time_frequency_measure( load_mats ...
  , 'load_func', @load ...
  , 'get_data_func', @(x) x.coh_mat ...
  , 'get_labels_func', @(x) x.coh_labs ...
  , 'get_time_func', @(x) x.coh_t ...
  , 'get_freqs_func', @(x) x.coh_f ...
  , 'transform_func', @(coh, labels, f, t) transform_per_trial_coh_do_norm_contrast(coh, labels, f, t, norm_across) ...
);

%%  social gaze coh

coh_p = fullfile( data_root, 'iti_gaze_aligned_coh' );
coh_mats = shared_utils.io.findmat( coh_p );

norm_each = { 'channel', 'regions', 'trialtypes', 'administration', 'looks_to' };
norm_across = { 'magnitudes', 'contexts', 'sessions', 'blocks', 'trials' };
tform_coh = @(coh, labels, f, t) transform_per_trial_coh_do_norm_contrast(coh, labels, f, t, norm_each, norm_across, false);

[coh, coh_labels, f, t] = bfw.load_time_frequency_measure( coh_mats ...
  , 'load_func', @shared_utils.io.fload ...
  , 'get_data_func', @(x) x.data ...
  , 'get_labels_func', @(x) x.labels ...
  , 'get_time_func', @(x) x.t ...
  , 'get_freqs_func', @(x) x.f ...
  , 'transform_func',  tform_coh ...
);

%%

[band_coh, band_coh_labs] = dsp3.get_band_means( ...
  coh, coh_labels', f, {[10, 20], [35, 51]}, {'beta', 'gamma'} );

%%
sub_each = setdiff( getcats(band_coh_labs), 'outcomes' );
[sn, sn_labs] = dsp3.summary_binary_op( ...
  band_coh, band_coh_labs', sub_each, 'self', 'none', @minus, @nanmean );
setcat( sn_labs, 'outcomes', 'self-none' );
[on, on_labs] = dsp3.summary_binary_op( ...
  band_coh, band_coh_labs', sub_each, 'other', 'none', @minus, @nanmean );
setcat( on_labs, 'outcomes', 'other-none' );

band_contrast_labs = [ sn_labs'; on_labs ];
band_contrast_coh = [ sn; on ];

%%  stats compare time window

plt_data = band_coh;
plt_labs = band_coh_labs';

plt_mask = pipe( rowmask(plt_labs) ...
  , @(m) find(plt_labs, {'acc_bla', 'pre', 'choice'}, m) ...
  , @(m) find(plt_labs, {'monkey'}, m) ...
  , @(m) find(plt_labs, {'beta', 'gamma'}, m) ...
);

t_win = [0.05, 0.35];
t_ind = t >= t_win(1) & t <= t_win(2);
plt_data = nanmean( plt_data(:, t_ind), 2 );

rs_each = { 'regions', 'looks_to', 'bands' };

anova_outs = dsp3.anova1( plt_data, plt_labs', rs_each, 'outcomes', 'mask', plt_mask );
ttest_outs = dsp3.ttest2( plt_data, plt_labs', rs_each, 'self-none', 'other-none', 'mask', plt_mask, 'descriptive_funcs', dsp3.nandescriptive_funcs );

sn_rs_outs = dsp3.ranksum( plt_data, plt_labs', rs_each, 'self', 'none', 'mask', plt_mask );
on_rs_outs = dsp3.ranksum( plt_data, plt_labs', rs_each, 'other', 'none', 'mask', plt_mask );

%%  

f_ind = f >= 10 & f <= 60;

over_freq_coh = nanmean( coh(:, f_ind, t >= 0.05 & t <= 0.35), 3 );
plt_labs = coh_labels';

plt_mask = pipe( rowmask(plt_labs) ...
  , @(m) find(plt_labs, {'acc_bla', 'pre', 'choice'}, m) ...
  , @(m) find(plt_labs, {'monkey'}, m) ...
);

pl = plotlabeled.make_common();
pl.add_smoothing = true;
pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.1);
pl.one_legend = false;

pl.x = f(f_ind);

[axs, hs, inds] = pl.lines( over_freq_coh(plt_mask, :), prune(plt_labs(plt_mask)) ...
  , 'outcomes', {'regions', 'looks_to', 'trialtypes', 'administration'} );

xlabel( axs(1), 'Frequency (Hz)' );
shared_utils.plot.set_ylims( axs, [-1.5e-2, 1.5e-2] );

bi = shared_utils.vector.slidebin( 1:size(over_freq_coh, 2), 3, 3, true );
is_sig = false( 1, size(over_freq_coh, 2) );

for i = 1:numel(bi)
  freq_coh = nanmean( over_freq_coh(plt_mask, bi{i}), 2 );  
  outs = dsp3.anova1( freq_coh, prune(plt_labs(plt_mask)), {}, 'outcomes' );
  is_sig(i) = outs.anova_tables{1}.Prob_F{1} < 0.05;
end

% is_sig = false( 1, size(over_freq_coh, 2) );
% for i = 1:size(over_freq_coh, 2)  
%   outs = dsp3.anova1( over_freq_coh(plt_mask, i), prune(plt_labs(plt_mask)), {}, 'outcomes' );
%   is_sig(i) = outs.anova_tables{1}.Prob_F{1} < 0.05;
% end

if ( sum(is_sig) > 0 )
  for i = 1:numel(axs)
    hold( axs(i), 'on' );
    x = f(f_ind);
    x = x(is_sig); 
    plot( axs(i), x, max(get(axs(i), 'ylim')) - diff(get(axs(i), 'ylim')) * 0.25, 'k*' );
  end
end

%%

plt_data = band_coh;
plt_labs = band_coh_labs';
% plt_mask = pipe( rowmask(plt_labs) ...
%   , @(m) find(plt_labs, {'spk-acc', 'lfp-bla', 'pre', 'choice'}, m) ...
%   , @(m) find(plt_labs, {'looks-to-monkey', 'looks-to-bottle'}, m) ...
%   , @(m) find(plt_labs, 'beta', m) ...
% );

plt_mask = pipe( rowmask(plt_labs) ...
  , @(m) find(plt_labs, {'acc_bla', 'pre', 'choice'}, m) ...
  , @(m) find(plt_labs, {'monkey'}, m) ...
  , @(m) find(plt_labs, {'beta', 'gamma'}, m) ...
);

t_ind = t >= -.150 & t <= .600;

plt_data = plt_data(plt_mask, t_ind);
plt_labs = prune( plt_labs(plt_mask) );

both_ind = find( plt_labs, 'looks-to-bottle-and-monkey' );
setcat( plt_labs, 'looks_to', 'looks-to-bottle', both_ind );
setcat( plt_labs, 'looks_to', 'looks-to-anywhere-else', find(plt_labs, 'looks-to-bottle') );

pl = plotlabeled.make_common();
pl.add_smoothing = true;
pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.3);
pl.one_legend = false;

pl.x = t(t_ind);

[axs, hs, inds] = pl.lines( plt_data, plt_labs, 'outcomes' ...
  , {'regions', 'bands', 'looks_to', 'trialtypes', 'administration'} );

% shared_utils.plot.set_ylims( axs, [-1.5e-2, 1.5e-2] );

% dsp3.compare_series( axs, inds, plt_data, @ranksum ...
%   , 'x', t(t_ind), 'series_handles', hs, 'fig', gcf );

figure(1);

%%  ff granger

granger_mats = shared_utils.io.findmat( '/Volumes/external3/data/changlab/ptp-vicarious-reward/f2fgrang-smaller' );

[granger, granger_labels, f, t] = bfw.load_time_frequency_measure( granger_mats ...
  , 'load_func', @load ...
  , 'get_data_func', @(x) replace_inf_with_nan(real(x.granger)) ...
  , 'get_labels_func', @(x) x.labels ...
  , 'get_time_func', @(x) x.t ...
  , 'get_freqs_func', @(x) x.f ...
);

%%

[band_granger, band_granger_labs] = dsp3.get_band_means( ...
  granger, granger_labels', f, {[10, 20], [35, 51]}, {'beta', 'gamma'} );

plt_data = band_granger;
plt_labs = band_granger_labs';

plt_mask = pipe( rowmask(plt_labs) ...
  , @(m) find(plt_labs, {'pre', 'self', 'none', 'other'}, m) ...
  , @(m) find(plt_labs, {'beta', 'gamma'}, m) ...
);

t_ind = true( size(t) );

plt_data = plt_data(plt_mask, t_ind);
plt_labs = prune( plt_labs(plt_mask) );

pl = plotlabeled.make_common();
pl.add_smoothing = true;
pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.5);
pl.one_legend = false;

pl.x = t(t_ind);

[axs, hs, inds] = pl.lines( plt_data, plt_labs, 'outcomes' ...
  , {'region', 'bands', 'trialtypes', 'administration'} );

% dsp3.compare_series( axs, inds, plt_data, @ranksum ...
%   , 'x', t(t_ind), 'series_handles', hs, 'fig', gcf );

figure(1);


%%  split by reward size

tot_mag_size_data = {};
tot_mag_size_labels = {};

for i = 1:size(sub_fs,1)
  
fprintf( '\n%d of %d', i, size(sub_fs, 1) );

loadStruct = load(sub_fs{i});

if strcmp(loadStruct.spkReg, loadStruct.lfpRegs)
  continue;
end

lfp_labs = loadStruct.ccLFPInfo;
lfp_labs(:, 2) = cellfun( @(x) sprintf('lfp-%s', x), lfp_labs(:, 2), 'un', 0 );
lfp_cols = { 'days', 'lfp-channels', 'magnitudes' ...
  , 'outcomes', 'trialtypes', 'administration', 'lfp-regions' };
lfp_ls = fcat.from( lfp_labs, lfp_cols );
assert( rows(lfp_ls) == size(loadStruct.cohMatrix, 3) );
spk_labs = loadStruct.ccSpikeInfo;
spk_labs(:, 3) = cellfun( @(x) sprintf('spk-%s', x), spk_labs(:, 2), 'un', 0 );
spk_ls = fcat.from( spk_labs, {'days', 'spk-channels', 'spk-regions', 'unit_index'} );
join( lfp_ls, spk_ls );

mu_each = { 'trialtypes', 'outcomes', 'administration', 'magnitudes' };
[mu_labs, mu_I, mu_C] = keepeach( lfp_ls', mu_each );

coh_mat = permute( loadStruct.cohMatrix, [3, 1, 2] );
mu_dat = bfw.row_nanmean( coh_mat, mu_I );

[ctx_labs, ctx_I, ctx_C] = keepeach( lfp_ls', {'trialtypes', 'administration', 'magnitudes'} );

visited = false( rows(mu_labs), 1 );

mag_size_labs = fcat();
mag_size_dat = [];
for j = 1:numel(ctx_I)
  sn_ind = find( lfp_ls, [ctx_C(1:2, j)', {'self', 'none'}] );
  on_ind = find( lfp_ls, [ctx_C(1:2, j)', {'other', 'none'}] );
  mu_sn = nanmean( coh_mat(sn_ind, :, :), 1 );
  mu_on = nanmean( coh_mat(on_ind, :, :), 1 );
  
  mu_s_ind = find( mu_labs, [ctx_C(:, j)', {'self'}] );
  mu_n_ind = find( mu_labs, [ctx_C(:, j)', {'none'}] );
  mu_o_ind = find( mu_labs, [ctx_C(:, j)', {'other'}] );
%   assert( (isempty(mu_s_ind) && isempty(mu_n_ind)) || ...
%     (numel(mu_s_ind) == 1 && numel(mu_n_ind) == 1) );
  
  mu_diff_sn = (mu_dat(mu_s_ind, :, :) - mu_sn) - (mu_dat(mu_n_ind, :, :) - mu_sn);
  mu_diff_on = (mu_dat(mu_o_ind, :, :) - mu_on) - (mu_dat(mu_n_ind, :, :) - mu_on);
  
  if ( isempty(mu_diff_sn) )
    mu_diff_sn = nan( 1, size(mu_dat, 2), size(mu_dat, 3) );
  end
  
  sb_labs = setcat( ctx_labs(j), 'outcomes', 'self-none' );
  append( mag_size_labs, sb_labs );
  mag_size_dat = [ mag_size_dat; mu_diff_sn ];
  
  if ( isempty(mu_diff_on) )
    mu_diff_on = nan( 1, size(mu_dat, 2), size(mu_dat, 3) );
  end
  
  on_labs = setcat( ctx_labs(j), 'outcomes', 'other-none' );
  append( mag_size_labs, on_labs );
  mag_size_dat = [ mag_size_dat; mu_diff_on ];
end

keep_f = loadStruct.coh_f <= 100;
assert( numel(keep_f) == size(mag_size_dat, 2) );

binned_t = shared_utils.vector.slidebin( loadStruct.t, 150, 50, true );
binned_t = cellfun( @median, binned_t );
assert( numel(binned_t) == size(mag_size_dat, 3) );
mag_size_dat = mag_size_dat(:, keep_f, :);

binned_f = loadStruct.coh_f(keep_f);
assert( numel(binned_f) == size(mag_size_dat, 2) );

assert_ispair( mag_size_dat, mag_size_labs );
tot_mag_size_labels{end+1} = mag_size_labs;
tot_mag_size_data{end+1} = mag_size_dat;

end

%%

mag_size_dat = vertcat( tot_mag_size_data{:} );
mag_size_labs = vertcat( fcat, tot_mag_size_labels{:} );
assert_ispair( mag_size_dat, mag_size_labs );

%%  

stored_coh = load( '/Volumes/external3/data/changlab/ptp-vicarious-reward/site-mean-sfcoh-by-magnitude/coh.mat' );
mag_size_dat = stored_coh.mag_size_dat;
mag_size_labs = stored_coh.mag_size_labs;
assert_ispair( mag_size_dat, mag_size_labs );
binned_t = stored_coh.binned_t;
binned_f = stored_coh.binned_f;

%%
f_win = [35, 51];
t_win = [50, 350];

% f_win = [10, 20];

t_ind = binned_t > t_win(1) & binned_t <= t_win(2);
f_ind = binned_f >= f_win(1) & binned_f <= f_win(2);
mu_dat = nanmean( nanmean(mag_size_dat(:, f_ind, t_ind), 2), 3 );

plt_mask = pipe( rowmask(mag_size_labs) ...
  , @(m) find(mag_size_labs, {'pre'}, m) ...
  , @(m) findnone(mag_size_labs, {'no_reward'}, m) ...
  , @(m) find(mag_size_labs, 'bla', m) ...
);

plt_subset = mu_dat(plt_mask, :, :);
plt_labs = prune( mag_size_labs(plt_mask) );

% each = setdiff( getcats(plt_labs), 'outcomes' );
% [plt_subset, plt_labs] = dsp3.summary_binary_op( ...
%   plt_subset, plt_labs', each, 'self-none', 'other-none', @minus, @nanmean );

comp_each = { 'trialtypes', 'administration', 'outcomes' };
[comp_labs, comp_I] = retaineach( plt_labs, comp_each );
addcat( comp_labs, 'comparison' );

clabs = fcat();
cps = [];
for i = 1:numel(comp_I)
  ci = comp_I{i};
  lowi = find( plt_labs, 'low', ci );
  medi = find( plt_labs, 'medium', ci );
  highi = find( plt_labs, 'high', ci );
  p_low_med = ranksum( plt_subset(lowi), plt_subset(medi) );
  p_med_high = ranksum( plt_subset(medi), plt_subset(highi) );
  append( clabs, comp_labs, i );
  addsetcat( clabs, 'comparison', 'low-v-med', rows(clabs) );
  append( clabs, comp_labs, i );
  addsetcat( clabs, 'comparison', 'med-v-high', rows(clabs) );
  cps = [ cps; p_low_med; p_med_high ];
end

comp_each = { 'trialtypes', 'administration', 'magnitudes' };
rs_outs = dsp3.ranksum( plt_subset, plt_labs, comp_each, 'self-none', 'other-none' );

%

pl = plotlabeled.make_common();
pl.x_order = { 'low', 'medium', 'high' };
axs = pl.errorbar( plt_subset, plt_labs, 'magnitudes', {'outcomes', 'trialtypes'}, {} );
set( axs, 'ylim', [-0.025, 0.025] );

figure( 1 );

%%

function lfp_ls = make_sf_labels(loadStruct)

lfp_labs = loadStruct.ccLFPInfo;
lfp_labs(:, end) = cellfun( @(x) sprintf('lfp-%s', x), lfp_labs(:, end), 'un', 0 );
lfp_cols = { 'days', 'lfp-channels', 'magnitudes' ...
  , 'outcomes', 'trialtypes', 'administration', 'lfp-regions' };
lfp_ls = fcat.from( lfp_labs, lfp_cols );
assert( rows(lfp_ls) == size(loadStruct.cohMatrix, 3) );
spk_labs = loadStruct.ccSpikeInfo;
spk_labs(:, 3) = cellfun( @(x) sprintf('spk-%s', x), spk_labs(:, 3), 'un', 0 );
spk_ls = fcat.from( spk_labs, {'days', 'spk-channels', 'spk-regions', 'unit_index'} );
join( lfp_ls, spk_ls );

end

function x = replace_inf_with_nan(x)

x(~isfinite(x)) = nan;

end
