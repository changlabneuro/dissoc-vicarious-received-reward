clearvars
clc
%% Paths to COH matrices
% path2CohMat = '/Users/putnampt/Dropbox (ChangLab)/ptp/VicariousRewardLFP/COH_MAT/COH_MATRIX_STDNORM.mat';
% path2CohMat = '/Volumes/external3/ptp-vicarious-reward/COH_MAT/COH_MATRIX_STDNORM.mat';
path2CohMat = '/Volumes/external3/data/changlab/ptp-vicarious-reward/COH_MAT/COH_MATRIX_STDNORM.mat';

% path2CohMat = '~/source/changlab/ptp-vicarious-reward-data/COH_MAT/COH_MATRIX_STDNORM.mat';

%% Settings
cnorm = true;
alphaVal  = 0.05;
%% Load SFC data
load(path2CohMat);
%% Post-load settings
minPlotFreq = 10;
maxPlotFreq = 60;
minPlotMs = -250;
maxPlotMs = 750;
minDisplayMs = -150;
maxDisplayMs = 600;
smoothFac = 3;
ROI_MinFreq = 35;
ROI_MaxFreq = 51;

% ROI_MinFreq = 10;
% ROI_MaxFreq = 20;

ROI_MinMs = 50;
ROI_MaxMs = 400;
saveFigs = true;

ttest_bin_by_bin = true;

if ( ttest_bin_by_bin )
  bin_p_func = @ttest_p;
  comp_p_func = @ttest2_p;
else
  bin_p_func = @signrank;
  comp_p_func = @ranksum;
end

%% Reduce size of averageMat
% averageMat: (frequencies x time x epoch x (acc->bla outcomes, bla -> acc outcomes) x trials

validFreqIdxs = find(loadStruct.coh_f >= minPlotFreq & loadStruct.coh_f < maxPlotFreq);
meanTimeBins = cellfun(@mean,loadStruct.binned_t)+loadStruct.t(1);
validTimeIdxs = find(meanTimeBins > minPlotMs & meanTimeBins < maxPlotMs);
averageMat = averageMat(validFreqIdxs, validTimeIdxs, :, :, :);
timeLabels = meanTimeBins(validTimeIdxs);
freqLabels = loadStruct.coh_f(validFreqIdxs);

%%  reformat data

[ptp_coh, ptp_labels] = linearize_average_mat( averageMat ...
  , shared_utils.io.fload(fullfile(fileparts(path2CohMat), 'site_info.mat')) ...
  , shared_utils.io.fload('trial_data.mat') ...
);

% [ptp_coh, ptp_labels] = linearize_average_mat( averageMat ...
%   , [] ...
%   , shared_utils.io.fload('trial_data.mat') ...
% );

%%  replot spectra

mask = pipe( rowmask(ptp_labels) ...
  , @(m) find(ptp_labels, {'acc_bla'}, m) ...
);

do_z_rescale = false;

tl = timeLabels;
fl = freqLabels;

t_ind = tl >= minDisplayMs & tl <= maxDisplayMs;
tl = tl(t_ind);

pl = plotlabeled.make_spectrogram( fl, tl );
pl.panel_order = {'choice'};
if ( do_z_rescale )
  pl.smooth_func = @z_rescale;
end

axs = pl.imagesc( ptp_coh(mask, :, t_ind), ptp_labels(mask) ...
  , {'direction', 'trialtype', 'outcome'} );

if ( do_z_rescale )
  shared_utils.plot.set_clims( axs, [-.014, 0.014] );
else
  shared_utils.plot.set_clims( axs, [-.01, 0.01] );
end

time_freq_axis_labels( axs, tl, fl );

for i = 1:numel(axs)
  time_freq_rectangle( axs(i), tl, fl, [50, 400], [35, 51] );
  time_freq_rectangle( axs(i), tl, fl, [50, 400], [10, 20] );
end
hold( axs, 'on' );
% shared_utils.plot.add_vertical_lines( axs, find(timeLabels==0) );

figure(1);

%%  replot lines

[band_means, band_labs] = dsp3.get_band_means(...
  ptp_coh, ptp_labels', freqLabels, {[10, 20], [35, 51]}, {'alpha_beta', 'gamma'}, @nanmean );

mask = pipe( rowmask(band_labs) ...
  , @(m) find(band_labs, {'bla_acc'}, m) ...
);

[I, id, C] = rowsets( 2, band_labs, {'trialtype', 'direction', 'bands'}, 'outcome', 'mask', mask );
[PI, PL] = plots.nest2( id, I, strrep(plots.cellstr_join(C), '_', ' ') );

clf;
axs = plots.panels( numel(PI), true );
tstats = cell( numel(axs), 1 );
for i = 1:numel(axs)
  axes( axs(i) );
  I = PI{i};
  colors = { 'r', 'b' };
  for j = 1:numel(I)
    data = band_means(I{j}, :);
    h = stdshade( data,.5, colors{j}, timeLabels, smoothFac );
    set( h, 'displayname', PL{i, 2}{j} );
    hold on;
  end
  if ( numel(I) == 2 )
    t_win = timeLabels > 50 & timeLabels < 400;
    mu0 = nanmean( band_means(I{1}, t_win), 2 );
    mu1 = nanmean( band_means(I{2}, t_win), 2 );
    [~, test_p, ~, tstats{i}] = ttest2( mu0, mu1 );
    [tstats{i}.mu0, tstats{i}.mu1] = deal( mu0, mu1 );
    if ( test_p < 0.05 )
      plot( axs(i), mean(timeLabels(t_win)), 1e-2, 'k*' );
    end
  end
  legend;
  xlim( axs(i), [minDisplayMs maxDisplayMs] );
  ylim( axs(i), [-1e-2, 1e-2] );
  title( PL{i, 1} );
  shared_utils.plot.add_horizontal_lines( axs(i), 0 );
end

figure( 1 );

%%  replot bars

[band_means, band_labs] = dsp3.get_band_means(...
  ptp_coh, ptp_labels', freqLabels, {[10, 20], [35, 51]}, {'alpha_beta', 'gamma'}, @nanmean );
band_means = nanmean( band_means(:, timeLabels > 50 & timeLabels < 400), 2 );

mask = pipe( rowmask(band_labs) ...
  , @(m) find(band_labs, 'choice', m) ...
);

figure(1);
[I, id, C] = rowsets( 3, band_labs, {'bands', 'outcome'}, 'direction', 'trialtype' ...
  , 'mask', mask );
[axs, ~, xs, PI] = plots.simplest_barsets( band_means, I, id ...
 , plots.strip_underscore(plots.cellstr_join(C)) ...
  , 'summary_func', @nanmean, 'error_func', @plotlabeled.nansem );

shared_utils.plot.match_ylims( axs );
ps = nan( numel(PI), 1 );
for i = 1:numel(PI)
  [~, p] = ttest2( band_means(PI{i}{1}), band_means(PI{i}{2}) );
  ps(i) = p;
  if ( p < 0.05 )
    x0 = mean( xs{i} );
    y = max(get(axs(i), 'ylim')) - diff(get(axs(i), 'ylim')) * 0.1;
    plot( axs(i), x0, y, 'k*' );
    text( axs(i), x0 + 0.1, y, sprintf('p = %0.4f', p) );
  end
end

%%  

dr = fileparts( fileparts(path2CohMat) );
mats = shared_utils.io.findmat( fullfile(dr, 'remade-sfcoh-contrast') );
[ptp_coh, ptp_labels, freqLabels, t] = bfw.load_time_frequency_measure( mats ...
  , 'get_data_func', @(x) x.coh ...
  , 'get_labels_func', @(x) x.labels ...
);

f_ind = freqLabels > minPlotFreq & freqLabels < maxPlotFreq;
ptp_coh = ptp_coh(:, f_ind, :);
freqLabels = freqLabels(f_ind);

renamecat( ptp_labels, 'trialtypes', 'trialtype' );
renamecat( ptp_labels, 'outcomes', 'outcome' );
direc = fcat.strjoin( ptp_labels(:, {'spk-region', 'lfp-region'}), 1 );
direc = strrep( direc, 'spk-acc_lfp-bla', 'acc_bla' );
direc = strrep( direc, 'spk-bla_lfp-acc', 'bla_acc' );
addsetcat( ptp_labels, 'direction', direc );
replace( ptp_labels, 'self-none', 'self-bottle' );
replace( ptp_labels, 'other-none', 'other-bottle' );

timeLabels = t * 1e3;

%%  time-freq specificity

mask = pipe( rowmask(ptp_labels) ...
  , @(m) find(ptp_labels, {'acc_bla', 'choice'}, m) ...
);

t_ind = timeLabels > 50 & timeLabels < 400;
t = timeLabels(t_ind);
t_coh = squeeze( nanmean(ptp_coh(:, :, t_ind), 3) );

pl = plotlabeled.make_common();
pl.x = freqLabels;
axs = pl.lines( ...
  t_coh(mask, :), ptp_labels(mask), 'outcome', {'trialtype', 'direction'} );

%%  time-freq spec

mask = pipe( rowmask(ptp_labels) ...
  , @(m) find(ptp_labels, {'acc_bla', 'choice'}, m) ...
);

l1 = 'self-bottle';
l2 = 'other-bottle';
f0 = 10;

% l1 = 'other-bottle';
% l2 = 'self-bottle';
% f0 = 35;

f_ind = find( freqLabels > f0 );

[test_labs, test_I] = retaineach( ptp_labels, {'trialtype', 'direction'}, mask );
out_I = eachcell( @(x) bfw.find_combinations(ptp_labels, {l1, l2}, x), test_I );
assert( numel(out_I) == 1 );

t_ind = find( timeLabels > 50 & timeLabels < 400 );
t_coh = squeeze( nanmean(ptp_coh(:, :, t_ind), 3) );

ps = nan( numel(out_I), numel(f_ind) );
for i = 1:numel(f_ind)
  f_mean = nanmean( t_coh(:, f_ind(1:i)), 2 );
  ps(:, i) = cellfun( @(x) select(2, @() ttest2(f_mean(x{1}), f_mean(x{2}), 'tail', 'right')), out_I );
end

pl = plotlabeled.make_common();
pl.x = freqLabels(f_ind);
axs = pl.lines( ps, test_labs, 'outcome', {'trialtype', 'direction'} );
hold( axs, 'on' );
ylim( axs, [0, 1] );
ylabel( axs(1), sprintf('p (%s > %s)', l1, l2) );
shared_utils.plot.add_horizontal_lines( axs, 0.05 );

%%  time-freq shift across time

mask = pipe( rowmask(ptp_labels) ...
  , @(m) find(ptp_labels, {'acc_bla', 'choice'}, m) ...
);

is_ab = true;
center_bin = true;

l1 = 'self-bottle';
l2 = 'other-bottle';
f_s = [10, 20];
t_s = [50, 400];

if ( ~is_ab )
  l1 = 'other-bottle';
  l2 = 'self-bottle';
  f_s = [35, 51];
end

f_win_size = sum( freqLabels > f_s(1) & freqLabels < f_s(2) );
t_win_size = sum( timeLabels > 50 & timeLabels < 400 );
% t_win_size = 5;

[test_labs, test_I] = retaineach( ptp_labels, {'trialtype', 'direction'}, mask );
out_I = eachcell( @(x) bfw.find_combinations(ptp_labels, {l1, l2}, x), test_I );
assert( numel(out_I) == 1 );

tot_ps = nan( numel(freqLabels), numel(timeLabels) );
for i = 1:numel(freqLabels)
  if ( 0 )
    f0 = max( 1, i - floor(f_win_size * 0.5) );  
  else  
    f0 = i;
  end
  f1 = min( numel(freqLabels), f0 + f_win_size - 1 );
  
  for j = 1:numel(timeLabels)
    if ( 0 )
      t0 = max( 1, j - floor(t_win_size * 0.5) );
    else
      t0 = j;
    end
    
    t1 = min( numel(timeLabels), t0 + t_win_size - 1 );
    base = squeeze( nanmean(ptp_coh(:, f0:f1, t0:t1), [2, 3]) );
    [~, p] = ttest2( base(out_I{1}{1}), base(out_I{1}{2}), 'tail', 'right' );
    tot_ps(i, j) = p;
  end
end

tLabel0 = timeLabels;
freqLabel0 = freqLabels;

if ( center_bin )
  tLabel0 = tLabel0 + uniquetol(diff(timeLabels)) * 0.5 * t_win_size;
  freqLabel0 = freqLabel0 + uniquetol(diff(freqLabels)) * 0.5 * f_win_size;
end

clf;
ax = gca;
imagesc( ax, tLabel0, freqLabel0, imgaussfilt(tot_ps, 0.75) );
[i, j] = find( tot_ps < 0.05 );
hold( ax, 'on' );
scatter( ax, tLabel0(j), freqLabel0(i), 'w*' );
ylabel( ax, 'freq' );
xlabel( ax, 'time (ms)' );
colorbar( ax );
title( sprintf('p (%s > %s)', l1, l2) );
shared_utils.plot.set_clims( ax, [0, 1] );

rectangle( ax, 'position', [t_s(1), f_s(1), diff(t_s), diff(f_s)] );

%%  time-freq spec across time

mask = pipe( rowmask(ptp_labels) ...
  , @(m) find(ptp_labels, {'acc_bla', 'choice'}, m) ...
);

l1 = 'self-bottle';
l2 = 'other-bottle';
f0 = 10;

% l1 = 'other-bottle';
% l2 = 'self-bottle';
% f0 = 35;

f_ind = find( freqLabels > f0 );

tot_ps = {};

[test_labs, test_I] = retaineach( ptp_labels, {'trialtype', 'direction'}, mask );
out_I = eachcell( @(x) bfw.find_combinations(ptp_labels, {l1, l2}, x), test_I );
assert( numel(out_I) == 1 );

t_ind = find( timeLabels > 50 );
for idx = 1:numel(t_ind)

t_coh = squeeze( nanmean(ptp_coh(:, :, t_ind(1:idx)), 3) );

ps = nan( numel(out_I), numel(f_ind) );
for i = 1:numel(f_ind)
  f_mean = nanmean( t_coh(:, f_ind(1:i)), 2 );
  ps(:, i) = cellfun( @(x) select(2, @() ttest2(f_mean(x{1}), f_mean(x{2}), 'tail', 'right')), out_I );
end

tot_ps{end+1, 1} = ps;
end

tot_ps = vertcat( tot_ps{:} );

clf;
ax = gca;
p_ind = tot_ps < 0.05;
imagesc( ax, timeLabels(t_ind), freqLabels(f_ind), tot_ps' );
[i, j] = find( tot_ps < 0.05 );
hold( ax, 'on' );
scatter( ax, timeLabels(t_ind(i)), freqLabels(f_ind(j)), 'w*' );
ylabel( ax, 'freq' );
xlabel( ax, 'time (ms)' );

%%  spectra per animal

mask = pipe( rowmask(ptp_labels) ...
  , @(m) find(ptp_labels, {'acc_bla', 'choice'}, m) ...
);

t_ind = true( size(timeLabels) );
t = timeLabels(t_ind);

hitch_ind = find( ptp_labels, {'hitch', 'self-bottle'}, mask );

mt_ind = timeLabels > 50 & timeLabels < 400;
f_ind = freqLabels > 35 & freqLabels < 51;
mean_coh = nanmean( ptp_coh(:, f_ind, mt_ind), [2, 3] );
above_thresh = mean_coh > 0.1 | mean_coh < -0.08;
within_thresh = mean_coh > -0.08;

% mask = intersect( hitch_ind, intersect(mask, find(within_thresh)) );

pl = plotlabeled.make_spectrogram( freqLabels, t );
axs = pl.imagesc( ptp_coh(mask, :, t_ind), ptp_labels(mask) ...
  , {'trialtype', 'outcome', 'direction', 'monkeys'} );
shared_utils.plot.set_clims( axs, [-2.0e-2, 2.0e-2] );
shared_utils.plot.fseries_yticks( axs, round(flip(freqLabels)), 2 );
shared_utils.plot.tseries_xticks( axs, round(timeLabels), 4 );

for i = 1:numel(axs)
  hold( axs(i), 'on' );
  time_freq_rectangle( axs(i), t, freqLabels, [50, 350], [10, 20] );
end

%%

mask = pipe( rowmask(ptp_labels) ...
  , @(m) find(ptp_labels, {'acc_bla', 'choice'}, m) ...
);

t_ind = timeLabels > 50 & timeLabels < 400;
f_ind = freqLabels > 10 & freqLabels < 20;

mean_coh = nanmean( ptp_coh(mask, f_ind, t_ind), [2, 3] );

pl = plotlabeled.make_common();
axs = pl.boxplot( mean_coh, ptp_labels(mask) ...
  , {'trialtype', 'outcome'}, {'trialtype', 'direction', 'monkeys'} );

% ylim( axs(1), [-0.5e-2, 0.5e-2] );

%%  decode one animal to the other

mask = pipe( rowmask(ptp_labels) ...
  , @(m) find(ptp_labels, {'acc_bla', 'choice'}, m) ...
);

[all_accs, dec_labs] = ptp_decode_by_animal( ptp_coh, ptp_labels, freqLabels, timeLabels, mask );

%%

pl = plotlabeled.make_common();
axs = pl.boxplot( all_accs, dec_labs ...
  , {'direction', 'trialtype'}, {'monkeys', 'band'} );
hold( axs, 'on' );
shared_utils.plot.add_horizontal_lines( axs, 0.5 );

%%

acc_I = findall( dec_labs, {'monkeys', 'band', 'trialtype'} );

%%

pl = plotlabeled.make_common();
pl.x = timeLabels(ti);
axs = pl.lines( all_accs, dec_labs ...
  , {'direction', 'trialtype'}, {'monkeys', 'band'} );
hold( axs, 'on' );
shared_utils.plot.add_horizontal_lines( axs, 0.5 );
shared_utils.plot.add_vertical_lines( axs, [50, 350] );

%%

mask = pipe( rowmask(ptp_labels) ...
  , @(m) find(ptp_labels, {'acc_bla', 'choice'}, m) ...
);

[all_accs, dec_labs] = ptp_decode_per_animal( ptp_coh, ptp_labels, freqLabels, timeLabels, mask );

%%  fraction of sites with consistent sign (other-bottle > self-bottle) per animal

mask = pipe( rowmask(ptp_labels) ...
  , @(m) find(ptp_labels, {'acc_bla', 'choice'}, m) ...
);

[sign_labels, I] = retaineach( ptp_labels, setdiff(getcats(ptp_labels), 'outcome'), mask );
sign_coh = nan( [numel(I), size(ptp_coh, [2, 3])] );
for i = 1:numel(I)
  sn_ind = find( ptp_labels, 'self-bottle', I{i} );
  on_ind = find( ptp_labels, 'other-bottle', I{i} );
  assert( isscalar(sn_ind) && isscalar(on_ind) );
  sign_coh(i, :, :) = sign( ptp_coh(on_ind, :, :) - ptp_coh(sn_ind, :, :) );
end

sign_coh(~isnan(sign_coh)) = sign_coh(~isnan(sign_coh)) > 0;

[p_labels, I] = retaineach( sign_labels, {'direction', 'monkeys', 'outcome', 'trialtype'} );
ps = nan( [numel(I), size(sign_coh, [2, 3])] );
for i = 1:numel(I)
  s = sum( sign_coh(I{i}, :, :), 1, 'omitnan' );
  d = numel(I{i}) - sum( isnan(sign_coh(I{i}, :, :)), 1 );
  ps(i, :, :) = s ./ d;
end

pl = plotlabeled.make_spectrogram( freqLabels, timeLabels );
axs = pl.imagesc( ps, p_labels, {'trialtype', 'outcome', 'direction', 'monkeys'} );

%%  time freq bins labeled by sign consistency

mask = pipe( rowmask(ptp_labels) ...
  , @(m) find(ptp_labels, {'acc_bla', 'choice'}, m) ...
);

[cons_coh, cons_labels, sign_coh, sign_labels] = ptp_sign_consistency( ptp_coh, ptp_labels, false, mask );

pl = plotlabeled.make_spectrogram( freqLabels, timeLabels );
pl.add_smoothing = false;
% pl.color_func = @(n) [[1, 0, 0]; [0, 0, 0]; [0, 0, 1]; [1, 1, 0]];
pl.color_func = @(n) [[1, 0, 0]; [0, 0, 0]; [0, 0, 1]; [1, 1, 0]];
axs = pl.imagesc( cons_coh, cons_labels, {'trialtype', 'outcome', 'direction', 'monkeys'} );

pl.color_func = @(n) [[0, 0, 0]; [1, 0, 0]];
axs = pl.imagesc( sign_coh+1, sign_labels, {'trialtype', 'outcome', 'direction', 'monkeys'} );

time_freq_axis_labels( axs, timeLabels, freqLabels );

for i = 1:numel(axs)
  hold( axs(i), 'on' );
  h = time_freq_rectangle( axs(i), timeLabels, freqLabels, [50, 350], [10, 20] );
  set( h, 'edgecolor', 'white' );
  set( h, 'linewidth', 4 );
end

figure(1);

%%  compare proportions of positive signs in tf window between animals

band_names = { 'alpha_beta', 'gamma' };
band_ranges = { [10, 20], [35, 51] };
% band_names = { 'gamma' };
% band_ranges = { [10, 20] };
band_labs = repset( addcat(sign_labels', 'band'), 'band', band_names );

mu_each = {'direction', 'monkeys', 'outcome', 'trialtype'};

props = [];
npos = [];
ns = [];
for i = 1:numel(band_names)
  for j = 1:size(sign_coh, 1)
    fsub = freqLabels > band_ranges{i}(1) & freqLabels < band_ranges{i}(2);
    tsub = timeLabels > 50 & timeLabels < 400;
    coh_sub = sign_coh(j, fsub, tsub);
    props(end+1, 1) = sum( coh_sub(:) ) / numel( coh_sub );
    ns(end+1, 1) = numel( coh_sub );
    npos(end+1, 1) = sum( coh_sub(:) );
  end
end

[prop_test_labels, I] = retaineach( ...
  band_labs, setdiff([{'band'}, mu_each], 'monkeys') );
ps = nan( numel(I), 1 );
for i = 1:numel(I)
  h_ind = find( band_labs, 'hitch', I{i} );
  k_ind = find( band_labs, 'kuro', I{i} );
  assert( isscalar(h_ind) && isscalar(k_ind) );
  [~, ps(i)] = prop_test( [npos(h_ind), npos(k_ind)], [ns(h_ind), ns(k_ind)], false );
end

%%  roi means

ab_f_ind = freqLabels > 10 & freqLabels < 24;
g_f_ind = freqLabels > 35 & freqLabels < 51;
t_ind = find( timeLabels > 50 & timeLabels < 400 );

% mean_dims = [2, 3];
mean_dims = [2];

ab_mean = squeeze( nanmean(ptp_coh(:, ab_f_ind, t_ind), mean_dims) );
g_mean = squeeze( nanmean(ptp_coh(:, g_f_ind, t_ind), mean_dims) );

combined_means = [ ab_mean; g_mean ];
combined_labs = repset( addcat(ptp_labels', 'bands'), 'bands', {'alpha_beta', 'gamma'} );

%%  site deltas

plt_mask = pipe( rowmask(combined_labs) ...
  , @(m) find(combined_labs, {'choice', 'acc_bla'}, m) ...
);

[site_labs, site_I] = retaineach( combined_labs, {'direction', 'site_index', 'bands'}, plt_mask );
site_deltas = nan( numel(site_I), size(combined_means, 2) );

for i = 1:numel(site_I)
  si = site_I{i};
  oni = find( combined_labs, 'other-bottle', si );
  sbi = find( combined_labs, 'self-bottle', si );
  assert( numel(oni) == numel(sbi) );
  is_miss = all( combined_means(oni, :) == 0 & combined_means(sbi, :) == 0, 1 );  
  site_deltas(i, :) = nanmean( combined_means(oni, :) - combined_means(sbi, :) );
  site_deltas(i, is_miss) = nan;
end

[I, C] = findeach( site_labs, {'monkeys', 'bands', 'direction', 'trialtype'} );
meds = zeros( numel(I), size(site_deltas, 2) );
for i = 1:numel(I)
  meds(i, :) = nanmean( site_deltas(I{i}, :), 1 );
end

%%  site bootstrapping, per animal

mask = pipe( rowmask(combined_labs) ...
  , @(m) find(combined_labs, {'choice', 'acc_bla'}, m) ...
);

[bs_labs, bs_I] = retaineach( combined_labs, {'direction', 'bands', 'monkeys'}, mask );
delta_sets = cell( size(bs_I) );
for i = 1:numel(bs_I)
  bsi = bs_I{i};
  s_I = findall( combined_labs, 'site_index', bsi );
  deltas = nan( 1e2, size(combined_means, 2) );
  for j = 1:1e2
    part = cvpartition( numel(s_I), 'holdout', 0.5 );
    grp = vertcat( s_I{part.training} );
    i0 = find( combined_labs, 'other-bottle', grp );
    i1 = find( combined_labs, 'self-bottle', grp );
    deltas(j, :) = nanmean( combined_means(i0, :) ) - nanmean( combined_means(i1, :) );
  end
  delta_sets{i} = deltas;
end

means = cate1( cellfun(@(x) mean(x, 1), delta_sets, 'un', 0) );
means = arrayfun( @(i) means(:, i), 1:size(means, 2), 'un', 0 );

t = table( bs_labs(:, 'monkeys'), means{:} ...
  , 'VariableNames', [{'monkeys'} ...
  , arrayfun(@(x) sprintf('t=%0.3f', timeLabels(t_ind(x))), 1:size(means, 2), 'un', 0)] );

%%

figure( 1 );
clf;

[I, id, C] = rowsets( 1, site_labs, {'monkeys', 'bands', 'direction', 'trialtype'} );
[PI, PL] = plots.nest1( id, I, plots.cellstr_join(C) );
axs = plots.panels( numel(PI) );

meds = nan( size(PI) );
tfs = false( size(PI) );
for i = 1:numel(PI)
  ax = axs(i);
  hist( ax, site_deltas(PI{i}{1}), 100 );
  title( ax, PL{i, 1} );
  [h, p] = ttest( site_deltas(PI{i}{1}) );
  tfs(i) = p < 0.05;
  meds(i) = nanmean( site_deltas(PI{i}{1}) );
end

shared_utils.plot.match_ylims( axs );

for i = 1:numel(axs)
  ld = diff( get(axs(i), 'ylim') );
  hold( axs(i), 'on' );
  text( axs(i), meds(i), max(get(axs(i), 'ylim')) - ld * 0.1 ...
    , sprintf('M = %0.3f', meds(i)) );
  shared_utils.plot.add_vertical_lines( axs(i), meds(i) );
end

%%  replot lines

plt_mask = pipe( rowmask(ptp_labels) ...
  , @(m) find(ptp_labels, {'choice', 'acc_bla'}, m) ...
);

freq_band = 'gamma';

switch ( freq_band )
  case 'alpha_beta'
    f_ind = freqLabels > 10 & freqLabels < 20;
  case 'gamma'
    f_ind = freqLabels > 35 & freqLabels < 51;
  otherwise
    error( 'Unrecognized freq band "%s".', freq_band );
end

f_mean = squeeze( nanmean(ptp_coh(:, f_ind, :), 2) );
addsetcat( ptp_labels, 'bands', freq_band );

figure( 1 );
pcats = { 'direction', 'trialtype', 'monkeys', 'bands' };
gcats = { 'outcome' };
[I, id, C] = rowsets( 2, ptp_labels, pcats, gcats, 'mask', plt_mask );

L = plots.strip_underscore( plots.cellstr_join(C) );
[PI, PL] = plots.nest2( id, I, L );
axs = plots.simplest_linesets( timeLabels, f_mean, PI, PL ...
  , 'summary_func', @(x) nanmean(x, 1) ...
  , 'error_func', @plotlabeled.nansem ...
  , 'smooth_func', @(x) cate1(...
      arrayfun(@(i) ptp_boxFilter(x(i, :), 3), 1:size(x, 1), 'un', 0)) ...
);

shared_utils.plot.match_ylims( axs );
shared_utils.plot.set_xlims( axs, [minDisplayMs, maxDisplayMs] );
shared_utils.plot.add_horizontal_lines( axs, 0 );

%%  lfp site frequencies
  
m_lfp = [SFCTable.lfp_day, SFCTable.lfp_chan, SFCTable.lfp_reg];
site_table_labs = fcat.from( categorical(m_lfp) ...
  , {'days', 'channels', 'regions'} );
cons = load( 'trial_data.mat' );
evt_labs = fcat.from( cons.consolidated.events.labels );
[~, m_C] = findall( evt_labs, {'days', 'monkeys'} );
match_I = bfw.find_combinations( site_table_labs, m_C(1, :) );

for i = 1:numel(match_I)
  if ( ~isempty(match_I{i}) )
    addsetcat( site_table_labs, 'monkeys', m_C(2, i), match_I{i} );
  end
end

site_table_labs = fcat.distinct( site_table_labs );
[ct_I, ct_C] = findall( site_table_labs, {'regions', 'monkeys'} );
freqs = cellfun( @numel, ct_I );

%%  cell frequencies

m_spk = [SFCTable.spk_day, SFCTable.spk_chan, SFCTable.spk_idx, SFCTable.spk_reg, SFCTable.spk_uuid];
site_table_labs = fcat.from( categorical(m_spk) ...
  , {'days', 'channels', 'unit_index', 'regions', 'unit_uuid'} );
site_table_labs = fcat.distinct( site_table_labs );

[~, m_C] = findall( evt_labs, {'days', 'monkeys'} );
match_I = bfw.find_combinations( site_table_labs, m_C(1, :) );

for i = 1:numel(match_I)
  if ( ~isempty(match_I{i}) )
    addsetcat( site_table_labs, 'monkeys', m_C(2, i), match_I{i} );
  end
end

[ct_I, ct_C] = findall( site_table_labs, {'regions', 'monkeys'} );
freqs = cellfun( @numel, ct_I );

%%  site pair frequencies

mat_subset = [ m_lfp, m_spk(:, 2:end) ];
mat_subset(:, 3) = cellfun( @(x) sprintf('lfp-%s', x), mat_subset(:, 3), 'un', 0 );
mat_subset(:, 6) = cellfun( @(x) sprintf('spk-%s', x), mat_subset(:, 6), 'un', 0 );

m_cols = { 'days', 'lfp_channel', 'lfp_region', 'spk_channel', 'unit_index', 'spk_region', 'unit_uuid' };
site_pair_labs = fcat.from( categorical(mat_subset), m_cols );

match_I = bfw.find_combinations( site_pair_labs, m_C(1, :) );

for i = 1:numel(match_I)
  if ( ~isempty(match_I{i}) )
    addsetcat( site_pair_labs, 'monkeys', m_C(2, i), match_I{i} );
  end
end

% acc_bla = find( site_pair_labs, {'spk-acc', 'lfp-bla'} );
acc_bla = find( site_pair_labs, {'lfp-acc', 'spk-bla'} );
hitch_ind = find( site_pair_labs, 'hitch', acc_bla );
kuro_ind = find( site_pair_labs, 'kuro', acc_bla );

%% Split apart into conditions and calculate mean

acc_bla = false;

t1_cued_normSelf_popCoh     = squeeze(averageMat(:,:,1,1,:));
t1_cued_normOther_popCoh    = squeeze(averageMat(:,:,1,2,:));
t1_choice_normSelf_popCoh    = squeeze(averageMat(:,:,1,3,:));
t1_choice_normOther_popCoh   = squeeze(averageMat(:,:,1,4,:));
t2_cued_normSelf_popCoh     = squeeze(averageMat(:,:,2,1,:));
t2_cued_normOther_popCoh    = squeeze(averageMat(:,:,2,2,:));
t2_choice_normSelf_popCoh    = squeeze(averageMat(:,:,2,3,:));
t2_choice_normOther_popCoh   = squeeze(averageMat(:,:,2,4,:));
t1_cued_normSelf_mCoh     = nanmean(t1_cued_normSelf_popCoh,3);
t1_cued_normOther_mCoh    = nanmean(t1_cued_normOther_popCoh,3);
t1_choice_normSelf_mCoh    = nanmean(t1_choice_normSelf_popCoh,3);
t1_choice_normOther_mCoh   = nanmean(t1_choice_normOther_popCoh,3);
t2_cued_normSelf_mCoh     = nanmean(t2_cued_normSelf_popCoh,3);
t2_cued_normOther_mCoh    = nanmean(t2_cued_normOther_popCoh,3);
t2_choice_normSelf_mCoh    = nanmean(t2_choice_normSelf_popCoh,3);
t2_choice_normOther_mCoh   = nanmean(t2_choice_normOther_popCoh,3);

if ( ~acc_bla )
  t1_choice_normSelf_mCoh = t2_choice_normSelf_mCoh;
  t1_choice_normOther_mCoh = t2_choice_normOther_mCoh;
  t1_choice_normSelf_popCoh = t2_choice_normSelf_popCoh;
  t1_choice_normOther_popCoh = t2_choice_normOther_popCoh;
  
  t1_cued_normSelf_mCoh = t2_cued_normSelf_mCoh;
  t1_cued_normOther_mCoh = t2_cued_normOther_mCoh;
  t1_cued_normSelf_popCoh = t2_cued_normSelf_popCoh;
  t1_cued_normOther_popCoh = t2_cued_normOther_popCoh;
end

%% Plot Figure 2 ACC --> BLA
fig2AH = figure(2);
clf
%%Choice
ROI_FreqIdxs = find(freqLabels > ROI_MinFreq & freqLabels < ROI_MaxFreq);
ROI_TimeIdxs= find(timeLabels > ROI_MinMs & timeLabels < ROI_MaxMs);
ROI_FreqRange =ROI_MaxFreq -ROI_MinFreq;
ROI_TimeRange = ROI_MaxMs-ROI_MinMs;
%% Choice Self SFC Spectro
axH(1) = subplot(2,4,1);
hold( axH(1), 'on' );
plotSpectroFuncZscaled(t1_choice_normSelf_mCoh, axH(1), timeLabels, freqLabels)
title('Self - Bottle');
xlabel('Time (ms)')
ylabel('Hz');
xlim([minDisplayMs maxDisplayMs])
rectangle('Position',...
  [ROI_MinMs ROI_MinFreq ROI_TimeRange ROI_FreqRange]...
  ,'EdgeColor','k',...
  'LineWidth',3,...
  'LineStyle', ':');
%% Choice Other SFC Spectro
axH(2) = subplot(2,4,2);
hold( axH(2), 'on' );
plotSpectroFuncZscaled(t1_choice_normOther_mCoh, axH(2), timeLabels, freqLabels)
title('Other - Bottle');
xlabel('Time (ms)')
ylabel('Hz');
rectangle('Position',...
  [ROI_MinMs ROI_MinFreq ROI_TimeRange ROI_FreqRange]...
  ,'EdgeColor','k',...
  'LineWidth',3,...
  'LineStyle', ':');
xlim([minDisplayMs maxDisplayMs])

%% Plotting SFC time course for choice trials
axH(3) = subplot(2,4,3:4);
cla( axH(3) );
hold on;
patch('XData', [ROI_MinMs ROI_TimeRange ROI_TimeRange ROI_MinMs],...
  'YData', [-1 -1 1 1],...
  'FaceColor','black',...
  'FaceAlpha',.1,...
  'LineStyle', 'none');
t1_choice_normSelf_ROICoh = t1_choice_normSelf_popCoh(ROI_FreqIdxs,:,:);
t1_choice_normSelf_ROICoh_avgXFrq = squeeze(nanmean(t1_choice_normSelf_popCoh(ROI_FreqIdxs,:,:),1));
t1_choice_normOther_ROICoh = t1_choice_normOther_popCoh(ROI_FreqIdxs,:,:);
t1_choice_normOther_ROICoh_avgXFrq = squeeze(nanmean(t1_choice_normOther_popCoh(ROI_FreqIdxs,:,:),1));
[choice_self_lineOut, ~] = stdshade(t1_choice_normSelf_ROICoh_avgXFrq',.5,'r', timeLabels, smoothFac);
[choice_other_lineOut, ~] = stdshade(t1_choice_normOther_ROICoh_avgXFrq',.5,'b', timeLabels, smoothFac);
axis tight
hline(0, 'k:');
vline(0, 'k-');
ntp = size(t1_choice_normSelf_ROICoh_avgXFrq,1);
pvals = nan(ntp,1);
diffs = nan( size(pvals) );

pvals_sb = nan( size(pvals) );
pvals_on = nan( size(pvals) );

for tp = 1:ntp
  tp_choice_normSelf = t1_choice_normSelf_ROICoh_avgXFrq(tp,:);
  tp_choice_normOther = t1_choice_normOther_ROICoh_avgXFrq(tp,:);
    
%   [pvals(tp),h,stats] = comp_p_func(tp_choice_normSelf,tp_choice_normOther);
  pvals(tp) = comp_p_func(tp_choice_normSelf,tp_choice_normOther);
  
  pvals_sb(tp) = bin_p_func( t1_choice_normSelf_ROICoh_avgXFrq(tp, :) );
  pvals_on(tp) = bin_p_func( t1_choice_normOther_ROICoh_avgXFrq(tp, :) );
end
[FDR] = mafdr(pvals,'BHFDR',true);
sigIdxs = find(pvals <= alphaVal);
% sigROIIdxs = intersect(ROI_TimeIdxs, sigIdxs);

sb_t_sub = nanmean( t1_choice_normSelf_ROICoh_avgXFrq(ROI_TimeIdxs, :), 1 );
ob_t_sub = nanmean( t1_choice_normOther_ROICoh_avgXFrq(ROI_TimeIdxs, :), 1 );
[~, p_twin] = ttest2( sb_t_sub(:), ob_t_sub(:) );

if ( 0 )
  
else
  sigY = (.003)*ones(size(sigIdxs));
  sigX = timeLabels(sigIdxs);
  scatter(sigX,sigY, 155, 'k*')
  
  sig_sb = find( pvals_sb <= alphaVal );
  sig_on = find( pvals_on <= alphaVal );
  
  sigY = (.008)*ones(size(sig_sb));
  sigX = timeLabels(sig_sb);
  scatter(sigX,sigY, 155, 'r*')
  
  sigY = (.005)*ones(size(sig_on));
  sigX = timeLabels(sig_on);
  scatter(sigX,sigY, 155, 'b*')
end


xlim([minDisplayMs maxDisplayMs])
ylim([-.01 0.01])
legend([choice_self_lineOut choice_other_lineOut],{'Self','Other'}, 'location', 'south')
title(sprintf('%3.1f to %3.1f Hz', ROI_MinFreq, ROI_MaxFreq));
xlabel('Time (ms)')
ylabel('Difference in Î³ coherence')


%%  compare time window

tp_choice_normSelf = nanmean( t1_choice_normSelf_ROICoh_avgXFrq(ROI_TimeIdxs, :), 1 );
tp_choice_normOther = nanmean( t1_choice_normOther_ROICoh_avgXFrq(ROI_TimeIdxs, :), 1 );
pval = ranksum( tp_choice_normSelf, tp_choice_normOther);

% p_sb = signrank( tp_choice_normSelf );
% p_on = signrank( tp_choice_normOther );
[~, p_sb] = ttest( tp_choice_normSelf );
[~, p_on] = ttest( tp_choice_normOther );

%%

plt = tp_choice_normSelf;
figure(1);
clf();
hist( plt, 40 );
mu = nanmean( plt );
hold( gca, 'on' );
shared_utils.plot.add_vertical_lines( gca, mu );
text( mu, max(get(gca, 'ylim')), sprintf('mu = %0.4f', mu) );

%% Cued Self SFC Spectro
axH(5) = subplot(2,4,5);
hold( axH(5), 'on' );
plotSpectroFuncZscaled(t1_cued_normSelf_mCoh, axH(5), timeLabels, freqLabels)
title('Self (Forced) - Bottle(Forced)');
xlabel('Time (ms)')
ylabel('Hz');
rectangle('Position',...
  [ROI_MinMs ROI_MinFreq ROI_TimeRange ROI_FreqRange]...
  ,'EdgeColor','k',...
  'LineWidth',3,...
  'LineStyle', ':');
xlim([minDisplayMs maxDisplayMs])
%% Cued Other SFC Spectro
axH(6) = subplot(2,4,6);
hold( axH(6), 'on' );
plotSpectroFuncZscaled(t1_cued_normOther_mCoh, axH(6), timeLabels, freqLabels)
title('Other(Forced) - Bottle(Forced)');
xlabel('Time (ms)')
ylabel('Hz');
rectangle('Position',...
  [ROI_MinMs ROI_MinFreq ROI_TimeRange ROI_FreqRange]...
  ,'EdgeColor','k',...
  'LineWidth',3,...
  'LineStyle', ':');
xlim([minDisplayMs maxDisplayMs])
%% Plotting SFC time course for forced trials
axH(8) = subplot(2,4,7:8);
cla( axH(8) );
hold on
patch('XData', [ROI_MinMs ROI_TimeRange ROI_TimeRange ROI_MinMs],...
  'YData', [-1 -1 1 1],...
  'FaceColor','black',...
  'FaceAlpha',.1,...
  'LineStyle', 'none');
t1_cued_normSelf_ROICoh = t1_cued_normSelf_popCoh(ROI_FreqIdxs,:,:);
t1_cued_normSelf_ROICoh_avgXFrq = squeeze(nanmean(t1_cued_normSelf_popCoh(ROI_FreqIdxs,:,:),1));
t1_cued_normOther_ROICoh = t1_cued_normOther_popCoh(ROI_FreqIdxs,:,:);
t1_cued_normOther_ROICoh_avgXFrq = squeeze(nanmean(t1_cued_normOther_popCoh(ROI_FreqIdxs,:,:),1));
[cued_self_lineOut, ~] = stdshade(t1_cued_normSelf_ROICoh_avgXFrq',.5,'r', timeLabels, smoothFac);
[cued_other_lineOut, ~] = stdshade(t1_cued_normOther_ROICoh_avgXFrq',.5,'b', timeLabels, smoothFac);
axis tight
hline(0, 'k:');
vline(0, 'k-');
xlabel('Time (ms)')
ylabel('Difference in Î³ coherence')
ntp = size(t1_cued_normSelf_ROICoh_avgXFrq,1);
pvals = nan(ntp,1);
pvals_sb = nan( size(pvals) );
pvals_on = nan( size(pvals_sb) );

for tp = 1:ntp
  tp_cued_normSelf = t1_cued_normSelf_ROICoh_avgXFrq(tp,:);
  tp_cued_normOther = t1_cued_normOther_ROICoh_avgXFrq(tp,:);
%   [pvals(tp),h,stats] = ranksum(tp_cued_normSelf,tp_cued_normOther);
  pvals(tp) = comp_p_func(tp_cued_normSelf,tp_cued_normOther);
  
  pvals_sb(tp) = bin_p_func( tp_cued_normSelf );
  pvals_on(tp) = bin_p_func( tp_cued_normOther );
end
[FDR] = mafdr(pvals,'BHFDR',true);
sigIdxs = find(FDR <= alphaVal);
sigROIIdxs = intersect(ROI_TimeIdxs, sigIdxs);

sb_t_sub = nanmean( t1_cued_normSelf_ROICoh_avgXFrq(ROI_TimeIdxs, :), 1 );
ob_t_sub = nanmean( t1_cued_normOther_ROICoh_avgXFrq(ROI_TimeIdxs, :), 1 );
[~, p_twin_cued] = ttest2( sb_t_sub(:), ob_t_sub(:) );

if ( 0 )
  sigROIIdxs = sigIdxs;
  sigY = (.008)*ones(size(sigROIIdxs));
  sigX = timeLabels(sigROIIdxs);
  scatter(sigX,sigY, 155, 'g*')
else
  sigY = (.003)*ones(size(sigIdxs));
  sigX = timeLabels(sigIdxs);
  scatter(sigX,sigY, 155, 'k*')
  
  sig_sb = find( pvals_sb <= alphaVal );
  sig_on = find( pvals_on <= alphaVal );
  
  sigY = (.008)*ones(size(sig_sb));
  sigX = timeLabels(sig_sb);
  scatter(sigX,sigY, 155, 'r*')
  
  sigY = (.005)*ones(size(sig_on));
  sigX = timeLabels(sig_on);
  scatter(sigX,sigY, 155, 'b*')
end


xlim([minDisplayMs maxDisplayMs])
ylim([-.01 0.01])
legend([cued_self_lineOut cued_other_lineOut],{'Self (Forced) - Bottle(Forced)','Other(Forced)-Bottle(Forced)'}, 'location', 'south')
title(sprintf('%3.1f to %3.1f Hz', ROI_MinFreq, ROI_MaxFreq));

%%

shared_utils.plot.set_clims( axH([1, 2, 5, 6]), [-.005, 0.005] );

%%

mu_diffs = nan( ntp, 1 );
mu_sn = nan( ntp, 1 );
mu_on = nan( ntp, 1 );
mu_p = nan( ntp, 1 );

for tp = 1:ntp
  tp_choice_normSelf = t1_choice_normSelf_ROICoh_avgXFrq(tp,:);
  tp_choice_normOther = t1_choice_normOther_ROICoh_avgXFrq(tp,:);
  mu_diffs(tp) = nanmean( tp_choice_normSelf ) - nanmean( tp_choice_normOther );
  
  mu_sn(tp) = nanmean( tp_choice_normSelf );
  mu_on(tp) = nanmean( tp_choice_normOther );
  mu_p(tp) = ttest2_p( tp_choice_normSelf, tp_choice_normOther );
end

figure(3);
plot( mu_sn, 'r' ); hold on;
plot( mu_on, 'b' );

sig_p = mu_p < 0.05;
scatter( find(sig_p), repmat(max(get(gca, 'ylim')), sum(sig_p), 1) );

%%

figure(3);
clf();
fl = freqLabels;
plot( timeLabels, nanmean(t1_cued_normSelf_mCoh(fl >= 35 & fl <= 50, :)) );
hold on;
plot( timeLabels, nanmean(t1_cued_normOther_mCoh(fl >= 35 & fl <= 50, :)) );
legend( {'self', 'other'} );
shared_utils.plot.add_horizontal_lines( gca, 0 );

%%

tp_cued_normSelf = nanmean( t1_cued_normSelf_ROICoh_avgXFrq(ROI_TimeIdxs, :), 1 );
tp_cued_normOther = nanmean( t1_cued_normOther_ROICoh_avgXFrq(ROI_TimeIdxs, :), 1 );
pval = ranksum( tp_cued_normSelf, tp_cued_normOther);

% p_sb = signrank( tp_choice_normSelf );
% p_on = signrank( tp_choice_normOther );
[~, p_sb] = ttest( tp_cued_normSelf );
[~, p_on] = ttest( tp_cued_normOther );

%%  bar plots for time window means

figure(1); clf;

sets = { {t1_choice_normSelf_ROICoh_avgXFrq, t1_choice_normOther_ROICoh_avgXFrq} ...
  , {t1_cued_normSelf_ROICoh_avgXFrq, t1_cued_normOther_ROICoh_avgXFrq} };

title_labs = { 'choice', 'cued' };
site_average_first = true;
include_lines = false;
means_of_smoothed = true;

for i = 1:numel(sets)

shp = ternary( include_lines, {2, 2}, {1, 2} );
bar_ax1 = subplot( shp{:}, i );

sb_time = sets{i}{1}(ROI_TimeIdxs, :);
ob_time = sets{i}{2}(ROI_TimeIdxs, :);

sb = nanmean( sb_time, 1 );
ob = nanmean( ob_time, 1 );
[~, p] = ttest2( sb, ob );

sem_sb = plotlabeled.nansem( sb' );
sem_ob = plotlabeled.nansem( ob' );

if ( site_average_first )
  if ( means_of_smoothed )
    sb_time = ptp_boxFilter( nanmean(sets{i}{1}, 2)', smoothFac )';
    ob_time = ptp_boxFilter( nanmean(sets{i}{2}, 2)', smoothFac )';
    store_sb_time = sb_time;
    store_ob_time = ob_time;
    
    sb_time = sb_time(ROI_TimeIdxs, :);
    ob_time = ob_time(ROI_TimeIdxs, :);
    bar( bar_ax1, [nanmean(sb_time, 1), nanmean(ob_time, 1)] );
  else
    bar( bar_ax1, [nanmean(sb_time, [2, 1]), nanmean(ob_time, [2, 1])] );
  end
else
  bar( bar_ax1, [nanmean(sb), nanmean(ob)] );
  hold( bar_ax1, 'on' );
  plot( bar_ax1, [1, 1], [nanmean(sb) - sem_sb * 0.5, nanmean(sb) + sem_sb * 0.5], 'k' );
  plot( bar_ax1, [2, 2], [nanmean(ob) - sem_ob * 0.5, nanmean(ob) + sem_ob * 0.5], 'k' );
end

set( bar_ax1, 'xticklabel', {'self-bottle', 'other-bottle'} );
set( bar_ax1, 'ylim', [-4.5e-3, 4.5e-3] );

if ( ~site_average_first )
  txt = ternary( p < 0.05, sprintf('p=%0.4f', p), 'ns' );
  text( bar_ax1, 1, max(get(bar_ax1, 'ylim')) - diff(get(bar_ax1, 'ylim')) * 0.1, txt );
end

title( bar_ax1, title_labs{i} );

if ( ~include_lines ), continue; end;

errs1 = plotlabeled.nansem( sets{i}{1}' );
errs2 = plotlabeled.nansem( sets{i}{2}' );
mus1 = nanmean(sets{i}{1}, 2)';
mus2 = nanmean(sets{i}{2}, 2)';

line_ax1 = subplot( 2, 2, i + 2 );
hold( line_ax1, 'on' );
h1 = plot( line_ax1, timeLabels, mus1, 'r' );
h2 = plot( line_ax1, timeLabels, mus2, 'b' );
plots.lineerrs( line_ax1, timeLabels, [mus1; mus2], [errs1; errs2], [[1, 0, 0];[0, 0, 1]] );
ylim( line_ax1, [-2e-2, 2e-2] );
shared_utils.plot.add_vertical_lines( line_ax1, [50, 350] );
shared_utils.plot.add_horizontal_lines( line_ax1, 0 );
legend( [h1, h2], {'self-bottle', 'other-bottle'} );

end


%% Figure 2 / Part B for inset
fig2BH = figure(3);
clf
axH(5) = subplot(2,4,4);
hold on;
choice_self_ROI_vals = squeeze(nanmean(nanmean(t1_choice_normSelf_popCoh(ROI_FreqIdxs, ROI_TimeIdxs, :),1),2));
choice_self_ROI_vals(isnan(choice_self_ROI_vals)) = [];
plotViolinPlot(choice_self_ROI_vals, 1, 0.01, [0 0 1], 'left')
choice_other_ROI_vals = squeeze(nanmean(nanmean(t1_choice_normOther_popCoh(ROI_FreqIdxs, ROI_TimeIdxs, :),1),2));
choice_other_ROI_vals(isnan(choice_other_ROI_vals)) = [];
plotViolinPlot(choice_other_ROI_vals, 1, 0.01, [1 0 0], 'right')
xticks([]);
xlim([.99 1.01])
ylim([-.1 .1])
p = ranksum(choice_self_ROI_vals,choice_other_ROI_vals);
title(sprintf('Wilcoxon rank sum p=%6.5f', p));
axH(9) = subplot(2,4,8);
hold on;
cued_self_ROI_vals = squeeze(nanmean(nanmean(t1_cued_normSelf_popCoh(ROI_FreqIdxs, ROI_TimeIdxs, :),1),2));
cued_self_ROI_vals(isnan(cued_self_ROI_vals)) = [];
plotViolinPlot(cued_self_ROI_vals, 1, 0.01, [0 0 1], 'left')
cued_other_ROI_vals = squeeze(nanmean(nanmean(t1_cued_normOther_popCoh(ROI_FreqIdxs, ROI_TimeIdxs, :),1),2));
cued_other_ROI_vals(isnan(cued_other_ROI_vals)) = [];
plotViolinPlot(cued_other_ROI_vals, 1, 0.01, [1 0 0], 'right')
xticks([]);
xlim([.99 1.01])
ylim([-.1 .1])
p = ranksum(cued_self_ROI_vals,cued_other_ROI_vals);
title(sprintf('Wilcoxon rank sum p=%6.5f', p));
%% Save Figure 2 / Part A
figure(2);
axes2Subset(1,1) = axH(1);
axes2Subset(2,1) = axH(2);
axes2Subset(3,1) = axH(5);
axes2Subset(4,1) = axH(6);
if cnorm
  [fig_2_cMin, fig_2_cMax] = normSubsetCAxes(axes2Subset);
end
manuallyRescaleCAxis(fig2AH, axes2Subset, fig_2_cMin, .023)
titleStr = sprintf('ACC Spikes-BLA FP - %3.1f to %3.1f Hz - %4d ms to %4d ms', ROI_MinFreq, ROI_MaxFreq, ROI_MinMs, ROI_MaxMs);
sgtitle(titleStr, 'FontSize', 14)
fig2AsaveStr = sprintf('Fig2_ACCSpikes-BLAFP_GammaSpectro_PartA');
fig2AsaveStr = [fig2AsaveStr, '.pdf'];
fig2ASavePath = fullfile(pwd, fig2AsaveStr);
makeLandscapePDF(fig2AH, fig2ASavePath);
%% Save Figure 2 / Part B
figure(3);
titleStr = sprintf('ACC Spikes-BLA FP - %3.1f to %3.1f Hz - %4d ms to %4d ms', ROI_MinFreq, ROI_MaxFreq, ROI_MinMs, ROI_MaxMs);
sgtitle(titleStr, 'FontSize', 14)
fig2BsaveStr = sprintf('Fig2_ACCSpikes-BLAFP_GammaSpectro_PartB');
fig2BsaveStr = [fig2BsaveStr, '.pdf'];
fig2BSavePath = fullfile(pwd, fig2BsaveStr);
makeLandscapePDF(fig2BH, fig2BSavePath);


%%

function p = ttest2_p(x, y)
[~, p] = ttest2( x, y );
end

function p = ttest_p(x)
[~, p] = ttest( x );
end

function y = z_rescale(x)

z_Mat = zscore(x);
f_Mat = imgaussfilt(z_Mat,2);
y = rescale(f_Mat,min(min(x)),max(max(x)));

end