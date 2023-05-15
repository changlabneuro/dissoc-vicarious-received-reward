dr = '/Volumes/external3/data/changlab/ptp-vicarious-reward';
% mats = shared_utils.io.findmat( fullfile(dr, 'high-res-sfcoh-contrast') );
mats = shared_utils.io.findmat( fullfile(dr, 'remade-sfcoh-contrast') );

%

[coh, labels, freqs, t] = bfw.load_time_frequency_measure( mats ...
  , 'get_data_func', @(x) x.coh ...
  , 'get_labels_func', @(x) x.labels ...
);

%%

site_each = {'lfp-channel', 'spk-channel', 'unit_id'};
mi = findall( labels, 'monkeys' );
n = cellfun( @(m) numel(findall(labels, site_each, m)), mi );
min_n = min( n );
bs_frac = 0.25;
iters = 1e1;

bs_coh_labs = fcat();
bs_coh = {};
for i = 1:iters  
  for j = 1:numel(mi)
    I = findall( labels, site_each, mi{j} );
    keep_s = sort( randperm(numel(I), max(1, floor(min_n * bs_frac))) );
    ind = vertcat( I{keep_s} );
    append( bs_coh_labs, labels, ind );
    bs_coh{end+1, 1} = coh(ind, :, :);
  end
end

bs_coh = vertcat( bs_coh{:} );

%%

m = find( labels, 'hitch' );
n = numel( findall(labels, {'lfp-channel', 'spk-channel', 'unit_id'}, m) );

%%

figure(1);

fs = [35, 51];
% fs = [10, 20];
ts = [0.05, 0.35];

plt_coh = coh;
plt_labels = labels;

plt_coh = bs_coh;
plt_labels = bs_coh_labs;

plt_mask = pipe( rowmask(plt_labels) ...
  , @(m) find(plt_labels, {'spk-acc', 'lfp-bla', 'choice'}, m) ...
);

f_ind = freqs >= 10 & freqs <= 60;
f = freqs(f_ind);
pl = plotlabeled.make_spectrogram( f, t );
axs = pl.imagesc( plt_coh(plt_mask, f_ind, :), plt_labels(plt_mask) ...
  , {'spk-region', 'lfp-region', 'outcomes', 'trialtypes', 'monkeys'} );

time_freq_axis_labels( axs, t*1e3, f );

hold( axs, 'on' );
shared_utils.plot.add_vertical_lines( axs, find(t==0) );
shared_utils.plot.set_clims( axs, [-1.5e-2, 1.5e-2] );

for i = 1:numel(axs)
  time_freq_rectangle( axs(i), t, f, [0.05, .350], [35, 51] );
  time_freq_rectangle( axs(i), t, f, [0.05, 0.35], [10, 20] );
end

%%

plt_mask = pipe( rowmask(coh) ...
  , @(m) find(labels, {'spk-acc', 'lfp-bla', 'choice'}, m) ...
);

bands = { [10, 20], [35, 51] };
band_names = { 'alpha_beta', 'gamma' };
coh_means = cate1( cellfun(...
  @(x) squeeze(nanmean(coh(plt_mask, freqs > x(1) & freqs < x(2), :), 2)), bands, 'un', 0) );

labs = repset( addcat(labels(plt_mask), 'band'), 'band', band_names );
assert_ispair( coh_means, labs );

pl = plotlabeled.make_common();
pl.add_smoothing = true;
pl.smooth_func = @(x) smoothdata( x, 'SmoothingFactor', 0.125 );
pl.x = t;
axs = pl.lines( coh_means, labs ...
  , {'outcomes'}, {'spk-region', 'lfp-region', 'trialtypes', 'monkeys', 'band'} );
hold( axs, 'on' );
shared_utils.plot.add_vertical_lines( axs, [0.05, 0.35] );
xlabel( axs, 'Time from reward (s)' );

%%

plt_mask = pipe( rowmask(coh) ...
  , @(m) find(labels, {'spk-acc', 'lfp-bla', 'choice'}, m) ...
);

t_ind = t >= 0.05 & t <= 0.35;
coh_mean = squeeze( nanmean(coh(:, :, t_ind), 3) );

pl = plotlabeled.make_common();
pl.add_smoothing = true;
pl.smooth_func = @(x) smoothdata(x, 'smoothingfactor', 0.5);
pl.x = freqs;
axs = pl.lines( coh_mean(plt_mask, :, :), labels(plt_mask) ...
  , {'outcomes'}, {'spk-region', 'lfp-region', 'trialtypes', 'monkeys'} );
hold( axs, 'on' );
shared_utils.plot.add_vertical_lines( axs, [10, 20, 35, 50] );
xlabel( axs, 'Hz' );