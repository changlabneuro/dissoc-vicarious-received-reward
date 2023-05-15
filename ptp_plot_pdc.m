dr = '/Volumes/external3/data/changlab/ptp-vicarious-reward';
% pdc_dir = fullfile( dr, 'pdc-orig-site-mean-mo-30' );
% pdc_dir = fullfile( dr, 'pdc-site-mean-mo-30' );
% pdc_dir = fullfile( dr, 'pdc-site-mean' );
% pdc_dir = fullfile( dr, 'dtf' );
pdc_dir = fullfile( dr, 'pdc-site-mean-mo-30' );
% pdc_dir = fullfile( dr, 'pdc-site-mean-nmatch-mo-30' );

[pdc, pdc_ls, pdc_f, pdc_t] = bfw.load_time_frequency_measure( ...
  shared_utils.io.findmat(pdc_dir) ...
  , 'load_func', @load ...
  , 'get_data_func', @(x) x.pdc ...
  , 'get_labels_func', @(x) x.pdc_ls ...
  , 'get_time_func', @(x) x.pdc_t ...
  , 'get_freqs_func', @(x) x.pdc_f ...
);

if ( 1 )
  % swap direction labels since they were coded incorrectly
  [dir_I, dir_C] = findall( pdc_ls, 'regions' );
  assert( numel(dir_I) == 2 );
  setcat( pdc_ls, 'regions', dir_C{2}, dir_I{1} );
  setcat( pdc_ls, 'regions', dir_C{1}, dir_I{2} );
end

%%

site_pairs = shared_utils.io.fload( fullfile(dr, 'site-pairs/pairs.mat') );
base_mask = find_site_pairs( site_pairs, pdc_ls );

outcome_contrast = false;

if ( 0 )
  debug_days = cellfun( @(x) sprintf('day__%s2017', x) ...
    , {'0104', '0105', '0106', '0107', '0108'}, 'un', 0 );
  base_mask = findor( pdc_ls, debug_days, base_mask );
end

pair_pdc = pdc(base_mask, :, :);
pair_pdc_ls = pdc_ls(base_mask);

if ( outcome_contrast )
  sub_each = {'days', 'sites', 'channels', 'regions', 'trialtypes', 'administration'};
  [sb, sblabs] = dsp3.summary_binary_op( pair_pdc, pair_pdc_ls', sub_each, 'self', 'none', @minus, @(x) nanmean(x, 1) );
  setcat( sblabs, 'outcomes', 'self-none' );
  [ob, oblabs] = dsp3.summary_binary_op( pair_pdc, pair_pdc_ls', sub_each, 'other', 'none', @minus, @(x) nanmean(x, 1) );  
  setcat( oblabs, 'outcomes', 'other-none' );
  pair_pdc = [ sb; ob ];
  pair_pdc_ls = [ sblabs'; oblabs ];
end

%%

if ( 0 )
  bands = pdc;
  band_labs = pdc_ls';
else
  bands = pair_pdc;
  band_labs = pair_pdc_ls';
end

band_indices = {pdc_f > 35 & pdc_f < 51, pdc_f > 10 & pdc_f < 20};
[bands, band_labs] = dsp3.get_band_means( bands, band_labs, pdc_f, band_indices, {'gamma', 'alpha_beta'} );

%%  lines over time

outs = ternary( outcome_contrast, {'self-none', 'other-none'}, {'self','none','other'} );
mask_fn = @(l, m) pipe( m ...
  , @(m) find(l, 'choice', m) ...
  , @(m) find(l, outs, m) ...
  , @(m) find(l, 'gamma', m) ...
);

mask = mask_fn( band_labs, rowmask(band_labs) );
plt_labs = band_labs(mask);
plt_pdc = bands(mask, :);

figure(2);
pl = plotlabeled.make_common();
pl.x = pdc_t;
pl.add_smoothing = true;
pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.6);
axs = pl.lines( plt_pdc, plt_labs, {'outcomes'}, {'trialtypes', 'regions', 'bands'} );

xlim( axs, [-500, 500] );
hold( axs, 'on' );
% ylim(axs, [0, 1.2]);
shared_utils.plot.add_vertical_lines( axs, [50, 350] );

%%  bar

assert_ispair( bands, band_labs );

t_win = pdc_t >= 50 & pdc_t <= 350;
tmean_pdc = squeeze( nanmean(bands(:, t_win), 2) );

mask_fn = @(l, m) pipe( m ...
  , @(m) find(l, 'choice', m) ...
  , @(m) find(l, outs, m) ...
  , @(m) find(l, 'pre', m) ...
);

mask = mask_fn( band_labs', rowmask(band_labs) );
plt_labs = band_labs(mask);
plt_pdc = tmean_pdc(mask, :);

figure(2);
pl = plotlabeled.make_common();
pl.x_order = {'self','other','none'};
[axs, inds] = pl.bar( plt_pdc, plt_labs, {'outcomes'}, 'regions', {'trialtypes', 'bands'} );

tps = cell( size(inds) );
for i = 1:numel(inds)
  tps{i} = nan( size(inds{i}, 1), 1 );
  for j = 1:size(inds{i}, 1)
    i0 = inds{i}{j, 1};
    i1 = inds{i}{j, 2};
    [~, tps{i}(j)] = ttest2( plt_pdc(i0), plt_pdc(i1) );
  end
end

%%  spectra

assert_ispair( pair_pdc, pair_pdc_ls );

mask_fn = @(l, m) pipe( m ...
  , @(m) find(l, 'choice', m) ...
  , @(m) find(l, outs, m) ...
  , @(m) find(l, 'pre', m) ...
);

sub_f = pdc_f < 100;

mask = mask_fn( pair_pdc_ls', rowmask(pair_pdc_ls) );
plt_labs = pair_pdc_ls(mask);
plt_pdc = pair_pdc(mask, sub_f, :);

figure(2);
pl = plotlabeled.make_spectrogram( pdc_f(sub_f), pdc_t );
pl.shape = [3, 2];
axs = pl.imagesc( plt_pdc, plt_labs, {'outcomes', 'trialtypes', 'regions'} );
shared_utils.plot.fseries_yticks( axs, flip(round(pdc_f(sub_f))), 5 );
shared_utils.plot.tseries_xticks( axs, pdc_t, 5 );

%%  lines over freq

assert_ispair( pair_pdc, pair_pdc_ls );

t_win = pdc_t >= 50 & pdc_t <= 350;
tmean_pdc = squeeze( nanmean(pair_pdc(:, :, t_win), 3) );

mask_fn = @(l, m) pipe( m ...
  , @(m) find(l, 'choice', m) ...
  , @(m) find(l, outs, m) ...
  , @(m) find(l, 'pre', m) ...
);

mask = mask_fn( pair_pdc_ls', rowmask(pair_pdc_ls) );
plt_labs = pair_pdc_ls(mask);
plt_pdc = tmean_pdc(mask, :);

figure(2);
pl = plotlabeled.make_common();
pl.x = pdc_f;
pl.add_smoothing = false;
pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.25);
axs = pl.lines( plt_pdc, plt_labs, {'outcomes'}, {'trialtypes', 'regions'} );
xlim( axs, [0, 100] );

hold( axs, 'on' );
shared_utils.plot.add_vertical_lines( axs, [35, 51] );

