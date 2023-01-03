%%

src_p = '/Volumes/external3/data/changlab/ptp-vicarious-reward/site-mean-sfcoh';
file_list = shared_utils.io.findmat( src_p );
[coh, coh_labels, f, t] = load_site_mean_coh( file_list );

%%

test_coh = load( '/Users/nick/Downloads/rwd_sfc_acc-bla_site_linearized.mat' );
coh = test_coh.mean_coh;
coh = coh(:, 1:numel(f), :);

coh_labels = test_coh.mean_labs';
renamecat( coh_labels, 'outcomes', 'outcome' );
renamecat( coh_labels, 'sites', 'site' );
addsetcat( coh_labels, 'administration', 'pre' );
addcat( coh_labels, {'trial-type', 'lfp-region', 'spk-region'} );

%%

sub_each = { 'site', 'administration', 'trial-type' };
[sb_coh, sb_labs] = dsp3.summary_binary_op( ...
  coh, coh_labels', sub_each, 'self', 'both', @minus, @(x) nanmean(x, 1) );
setcat( sb_labs, 'outcome', 'self-both' );
[on_coh, on_labs] = dsp3.summary_binary_op( ...
  coh, coh_labels', sub_each, 'other', 'none', @minus, @(x) nanmean(x, 1) );
setcat( on_labs, 'outcome', 'other-none' );

contrast_coh = [ sb_coh; on_coh ];
contrast_labs = [ sb_labs'; on_labs ];

%%  spectra

plt_dat = contrast_coh;
plt_labs = contrast_labs';
plt_t = t;
plt_f = f;

plt_t = plt_t(1:size(plt_dat, 3));

t_ind = plt_t >= -200 & plt_t <= 600;
% t_ind = true( size(plt_t) );

f_ind = plt_f <= 70;
plt_dat = plt_dat(:, f_ind, t_ind);
plt_t = plt_t(t_ind);
plt_f = plt_f(f_ind);

mask = pipe( rowmask(plt_labs) ...
  , @(m) find(plt_labs, 'pre', m) ...
  , @(m) findor(plt_labs, {'lfp-bla-spk-acc', 'lfp-acc-spk-bla'}, m) ...
);

plt_dat = plt_dat(mask, :, :);
plt_labs = plt_labs(mask);

pl = plotlabeled.make_spectrogram( plt_f, plt_t );
% pl.smooth_func = @(x) imgaussfilt(x, 1.5);

axs = pl.imagesc( plt_dat, plt_labs, {'outcome', 'trial-type', 'site-pair'} );
shared_utils.plot.set_clims( axs, [-2.5e-2, 1.75e-2] );
shared_utils.plot.tseries_xticks( axs, round(plt_t), 5 );
shared_utils.plot.fseries_yticks( axs, round(flip(plt_f)), 5 );

%%

plt_dat = contrast_coh;
plt_labs = contrast_labs';
plt_t = t;
plt_f = f;

mask = pipe( rowmask(plt_labs) ...
  , @(m) find(plt_labs, 'pre', m) ...
  , @(m) findor(plt_labs, {'lfp-bla-spk-acc', 'lfp-acc-spk-bla'}, m) ...
);

plt_dat = plt_dat(mask, :, :);
plt_labs = plt_labs(mask);

% f_win = [35, 50];
f_win = [10, 20];
f_ind = plt_f >= f_win(1) & plt_f <= f_win(2);
plt_dat = squeeze( nanmean(plt_dat(:, f_ind, :), 2) );

t_win = [ -150, 650 ];
t_ind = plt_t >= t_win(1) & plt_t <= t_win(2);
plt_dat = plt_dat(:, t_ind);
plt_t = plt_t(t_ind);

pl = plotlabeled.make_common();
pl.x = plt_t;
[axs, ~, inds] = pl.lines( plt_dat, plt_labs, {'outcome'}, {'site-pair', 'trial-type'} );
shared_utils.plot.set_ylims( axs, [-4e-2, 4e-2] );

dsp3.compare_series( axs, inds, plt_dat, @(x, y) ranksum(x, y), 'x', plt_t );

%%

function labels = add_site_pair(labels)

[reg_I, reg_C] = findall( labels, {'lfp-region', 'spk-region'} );
for i = 1:numel(reg_I)
  site_pair = sprintf( '%s-%s', reg_C{1, i}, reg_C{2, i} );
  addsetcat( labels, 'site-pair', site_pair, reg_I{i} );
end

end

function [coh, labels, f, t] = load_site_mean_coh(files)

labels = cell( size(files) );
coh = cell( size(files) );
f = [];
t = [];

for i = 1:numel(files)
  shared_utils.general.progress( i, numel(files) );
  
  coh_f = shared_utils.io.fload( files{i} );
  coh{i} = coh_f.coh;
  labels{i} = coh_f.labels;
  add_site_pair( labels{i} );
  
  if ( ~hascat(labels{i}, 'site') && ~isempty(labels{i}) )
    addsetcat( labels{i}, 'site', sprintf('site-%d', i) );
  end
  
  if ( i == 1 )
    f = coh_f.f;
    t = coh_f.t;
  end
end

labels = vertcat( fcat, labels{:} );
coh = vertcat( coh{:} );

end