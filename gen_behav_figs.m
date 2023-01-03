%%  figure 1.

consolidated_trials = load( 'trial_data.mat' );

events = consolidated_trials.consolidated.events;
event_labels = fcat.from( events.labels );
event_ts = events.data;
reward_ts = event_ts(:, consolidated_trials.consolidated.event_key('rwdOn'));
start_ts = event_ts(:, consolidated_trials.consolidated.event_key('fixOn'));

mask = find( event_labels, 'pre' );
mask = findnone( event_labels, 'errors', mask );

[day_labs, day_I, day_C] = keepeach( event_labels', {'days', 'monkeys'}, mask );
freqs = cellfun( @numel, day_I );
[mean_labs, mean_I] = keepeach( day_labs', 'monkeys' );
mean_freqs = rowifun( @nanmean, mean_I, freqs );
dev_freqs = rowifun( @nanstd, mean_I, freqs );

%%  preference index

[day_labs, day_I] = keepeach( event_labels', 'days' );
prefs = zeros( numel(day_I), 2 );

for i = 1:numel(day_I)
  di = day_I{i};
  n_self = double( count(event_labels, 'self', di) );
  n_both = double( count(event_labels, 'both', di) );
  n_other = double( count(event_labels, 'other', di) );
  n_none = double( count(event_labels, 'none', di) );
  bs = (n_both - n_self) ./ (n_both + n_self);
  on = (n_other - n_none) ./ (n_other + n_none);
  prefs(i, 1) = bs;
  prefs(i, 2) = on;
end

mean_prefs = mean( prefs, 1 );
std_prefs = std( prefs, [], 1 );

p_bs = signrank( prefs(:, 1) );
p_on = signrank( prefs(:, 2) );

%%

gd = load( 'fig1_gazeData.mat' );
partner_roi = gd.rois.monkey;

%%

int_p = '/Volumes/external3/data/changlab/dsp3/intermediates';
xs = shared_utils.io.findmat( fullfile(int_p, 'gaze/x') );
% xs = xs(1:2);

roi = gd.rois.monkey;

ibs = cell( numel(xs), 1 );
ib_labels = cell( size(ibs) );
lb_s = 400e-3;
la_s = 2600e-3;
ib_t_series = linspace( -lb_s, la_s, size(summary_dat, 2) );

parfor i = 1:numel(xs)
  fprintf( '\n %d of %d', i, numel(xs) );
  
  fname = shared_utils.io.filenames( xs{i}, true );
  xf = shared_utils.io.fload( xs{i} );
  yf = shared_utils.io.fload( fullfile(int_p, 'gaze/y', fname) );
  tf = shared_utils.io.fload( fullfile(int_p, 'gaze/t', fname) );
  ls = xf.labels;
  
  l_mask = rowmask( ls );
  e_mask = find( event_labels, combs(ls, 'days') );
  assert( numel(e_mask) == numel(l_mask) );
  
%   ls = event_labels(e_mask);
  e_times = reward_ts(e_mask);
  s_times = start_ts(e_mask);
  
  ib_x = xf.data >= partner_roi(1) & xf.data <= partner_roi(3);
  ib_y = yf.data >= partner_roi(2) & yf.data <= partner_roi(4);
  ib = ib_x & ib_y;
  
  t_series = tf.data(l_mask, :);
  t_series(:, 1) = 0;
  assert( size(t_series, 1) == numel(e_times) );
  
  fs = 2e-3;
  lb = lb_s / fs;
  la = la_s / fs;
  
  ib_mat = false( size(ib, 1), lb + la + 1 );
  for j = 1:numel(e_times)
    if ( e_times(j) == 0 ), continue; end    
    trow = t_series(j, :) + s_times(j);
    [~, mini] = min( abs(e_times(j) - trow) );
    ib_subset = ib(j, mini-lb:mini+la);
    ib_mat(l_mask(j), :) = ib_subset;
  end
  
  if ( 1 )
    keep_l = ls';
    mean_ibs = ib_mat;
  else
    [keep_l, summary_I] = keepeach( ...
      ls', {'days', 'outcomes', 'administration', 'trialtypes'}, l_mask );
    mean_ibs = rowifun( @(x) mean(x, 1), summary_I, ib_mat, 'un', 0 );
    mean_ibs = cate1( mean_ibs );
  end
  
  ibs{i} = mean_ibs;
  ib_labels{i} = keep_l';
end

ibs = vertcat( ibs{:} );
ib_labels = vertcat( fcat, ib_labels{:} );
assert_ispair( ibs, ib_labels );

%%

ib_mask = pipe( rowmask(ib_labels) ...
  , @(m) find(ib_labels, 'pre', m) ...
  , @(m) findnone(ib_labels, 'errors', m) ...
  , @(m) find(ib_labels, 'choice', m) ...
);

if ( 1 )
  summary_dat = double( ibs(ib_mask, :) );
  summary_labs = ib_labels(ib_mask);
else
  [summary_labs, summary_I] = keepeach( ib_labels', {'days', 'outcomes', 'trialtypes'}, ib_mask );
  summary_dat = cate1( rowifun(@(x) mean(x, 1), summary_I, ibs, 'un', 0) );
end

%%

pl = plotlabeled.make_common();
pl.smooth_func = @(s) smoothdata(s, 'movmean', 10);
pl.add_smoothing = true;
pl.group_order = { 'self', 'both', 'other' };
pl.x = ib_t_series;

[axs, hs, inds] = pl.lines( summary_dat, summary_labs, 'outcomes', 'trialtypes' );

for i = 1:numel(inds)
  ii = inds{i};
  line_labs = arrayfun( @(x) get(x, 'displayname'), hs{i}, 'un', 0 );
  s_ind = ii{contains(line_labs, 'self')};
  b_ind = ii{contains(line_labs, 'both')};
  o_ind = ii{contains( line_labs, 'other')};
  n_ind = ii{contains( line_labs, 'none')};
  
  sb_ps = ranksum_matrix( summary_dat, {s_ind}, {b_ind} );
  on_ps = ranksum_matrix( summary_dat, {o_ind}, {n_ind} );
  hold( axs(i), 'on' );
  
  lims = get( axs(i), 'ylim' );
  sb_sig = sb_ps < 0.05;  
  on_sig = on_ps < 0.05;
  h1 = scatter( axs(i), pl.x(sb_sig), repmat(max(lims), sum(sb_sig), 1), 'filled', 's' );
  h2 = scatter( axs(i), pl.x(on_sig), repmat(max(lims) - diff(lims) * 0.05, sum(on_sig), 1), 'filled', 'o' );
end

figure(1);

%%

[sem_day_labs, day_I] = keepeach( summary_labs', {'days', 'trialtypes'} );
day_sems = cate1( rowifun(@(x) mean(plotlabeled.sem(x)), day_I, summary_dat, 'un', 0) );
n_sems = 4;
sb_ts = nan( numel(day_I), 1 );
on_ts = nan( numel(day_I), 1 );

for i = 1:numel(day_I)
  di = day_I{i};
  
  s_ind = find( summary_labs, 'self', di );
  b_ind = find( summary_labs, 'both', di );
  o_ind = find( summary_labs, 'other', di );
  n_ind = find( summary_labs, 'none', di );
  
  sem_func = @(x) plotlabeled.sem(x);
  % sem_func = @(x) nanstd(x, [], 1);
  
  sb_day_sem = mean(sem_func(summary_dat([s_ind; b_ind], :)));
  on_day_sem = mean(sem_func(summary_dat([o_ind; n_ind], :)));
  
  sb_abs_mean_diffs = abs( ...
    mean(summary_dat(s_ind, :), 1) - mean(summary_dat(b_ind, :), 1) );
  on_abs_mean_diffs = abs( ...
    mean(summary_dat(o_ind, :), 1) - mean(summary_dat(n_ind, :), 1) );
  
  sb_above_thresh = find( sb_abs_mean_diffs > sb_day_sem * n_sems, 1 );
  on_above_thresh = find( on_abs_mean_diffs > on_day_sem * n_sems, 1 );
  
  if ( ~isempty(sb_above_thresh) )
    sb_ts(i) = ib_t_series(sb_above_thresh);
  end
  if ( ~isempty(on_above_thresh) )
    on_ts(i) = ib_t_series(on_above_thresh);
  end
end

[p, stats] = ranksum( sb_ts, on_ts );
nanmean( on_ts )
nanmean( sb_ts )

