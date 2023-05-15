data_root = '/Volumes/external3/data/changlab/ptp-vicarious-reward';

%%  load original normalized coherence matrix and split by new cell cluster labels

[coh, coh_labels, f, t] = load_new_split( data_root, false );

%%  load original spike field coherence matrix split by original cluster labels

[coh, coh_labels, f, t] = load_orig_split( data_root );

%%

smoothFac = 0;
timeLabels = t;
timeLabels = t - 25.5;

%%

% freq_win = [10, 20];
freq_win = [35, 51];

f_ind = f >= freq_win(1) & f < freq_win(2);
t_sub = timeLabels >= 50 & timeLabels <= 350;
f_mean = squeeze( nanmean(coh(:, f_ind, t_sub), [2, 3]) );

n = 1e2;
is_sig = false( n, 1 );

self_gt = true;
allow_replace = true;
up_sample = false;

for i = 1:n
  sampled = find_cell_cluster_subsets_matching_n( ...
    coh_labels, rowmask(coh_labels), allow_replace, up_sample );
%   sampled = rowmask( coh_labels );
  mask = pipe( sampled ...
    , @(m) find(coh_labels, 'acc_bla', m) ...
    , @(m) find(coh_labels, {'choice'}, m) ... 
    , @(m) find(coh_labels, {'cell-cluster2'}, m) ...
  );
  sampled(~ismember(sampled, mask)) = [];

  if ( self_gt )
    [al, bl] = deal( 'self-bottle', 'other-bottle' );
  else
    [al, bl] = deal( 'other-bottle', 'self-bottle' );
  end
  
  ai = find( coh_labels, al );
  bi = find( coh_labels, bl );
  ai(~ismember(ai, mask)) = [];
  bi(~ismember(bi, mask)) = [];
  
  real_mu_a = nanmean( f_mean(ai) );
  real_mu_b = nanmean( f_mean(bi) );
  
  num_ai = numel( ai );
  inds = [ ai; bi ];
  inds = inds(randperm(numel(inds)));
  null_ai = inds(1:num_ai);
  null_bi = inds(num_ai+1:end);
  
  null_mu_a = nanmean( f_mean(null_ai) );
  null_mu_b = nanmean( f_mean(null_bi) );
  
%   is_sig(i) = real_mu_a > real_mu_b;
%   is_sig(i) = real_mu_a > real_mu_b & ...
%     (null_mu_a < null_mu_b | (real_mu_a - real_mu_b > null_mu_a - null_mu_b));
  [~, p] = ttest2( f_mean(ai), f_mean(bi) );
  is_sig(i) = p < 0.05;
end

p_is_sig = 1 - sum( is_sig ) / numel( is_sig );

%%

figure(3);

% freq_win = [35, 51];
freq_win = [10, 20];
f_ind = f >= freq_win(1) & f < freq_win(2);
f_mean = squeeze( nanmean(coh(:, f_ind, :), 2) );

if ( 0 )
  base_mask = sampled;
else
  base_mask = rowmask( f_mean );
end

plt_labels = coh_labels';
plt_mask = pipe( base_mask ...
  , @(m) find(plt_labels, 'acc_bla', m) ...
  , @(m) find(plt_labels, {'choice'}, m) ... 
  , @(m) find(plt_labels, {'cell-cluster1', 'cell-cluster2'}, m) ...
);

[I, id, C] = rowsets( 2, plt_labels, {'direction', 'cell-cluster'}, {'outcome'}, 'mask', plt_mask );
[PI, PL] = plots.nest2( id, I, plots.cellstr_join(C) );
axs = plots.simple2( PI, PL, @(ax, inds, j, label) stdshade_cb(f_mean, timeLabels, inds, j, label) );
shared_utils.plot.add_horizontal_lines( axs, 0 );
shared_utils.plot.set_xlims( axs, [-150, 600] );

% shared_utils.plot.set_ylims( axs, [-1.e-2, 1.e-2] );
shared_utils.plot.set_ylims( axs, [-4.5e-2, 4.5e-2] );
shared_utils.plot.add_horizontal_lines( axs, [-1e-2, 1e-2] );
shared_utils.plot.add_vertical_lines( axs, [50, 350] );

lims = get( axs(1), 'ylim' );

t_sub = [ timeLabels >= 50 & timeLabels <= 350 ];

ps = nan( numel(PI), 1 );
for i = 1:numel(PI)
  ai = PI{i}{1};
  bi = PI{i}{2};
  bin_p = arrayfun( @(j) ttest2p(f_mean(ai, j), f_mean(bi, j)), 1:size(f_mean, 2) );
  sig_p = find( bin_p < 0.05 );
  scatter( axs(i), timeLabels(sig_p), repmat(max(lims), numel(sig_p), 1), 'k*' );
  
  [~, ps(i)] = ttest2( nanmean(f_mean(ai, t_sub), 2), nanmean(f_mean(bi, t_sub), 2) );
end

for i = 1:2
  col = ternary( i == 1, 'r*', 'b*' );
  ai = PI{1}{i};
  bi = PI{2}{i};
  
  bin_p = arrayfun( @(j) ttest2p(f_mean(ai, j), f_mean(bi, j)), 1:size(f_mean, 2) );
  sig_p = find( bin_p < 0.05 );
  for j = 1:numel(axs)
    scatter( axs(j), timeLabels(sig_p), repmat(max(lims) - diff(lims)*0.1*i, numel(sig_p), 1), col );
  end
end


site_mean_first = true;
means_of_smoothed = true;
figure(2);
clf();
axs = plots.panels( numel(PI), true );
for i = 1:numel(PI)
  hold( axs(i), 'on' );
  
  ai = PI{i}{1};
  bi = PI{i}{2};
  
  if ( site_mean_first )
    if ( means_of_smoothed )
      s0 = ptp_boxFilter( nanmean(f_mean(ai, :), 1), 3 );
      s1 = ptp_boxFilter( nanmean(f_mean(bi, :), 1), 3 );
      s0 = nanmean( s0(:, t_sub) );
      s1 = nanmean( s1(:, t_sub) );
    else
      s0 = nanmean( f_mean(ai, t_sub), [1, 2] );
      s1 = nanmean( f_mean(bi, t_sub), [1, 2] );
    end
    [~, x] = plots.bars( axs(i), [s0, s1]', PL{i, 2}, {''}, strrep(PL{i, 1}, '_', ' ') );
    set( axs(i), 'ylim', [-30e-3, 30e-3] );
  else
    s0 = nanmean( f_mean(ai, t_sub), 2 );
    s1 = nanmean( f_mean(bi, t_sub), 2 );
    [~, p] = ttest2( s0, s1 );

    p_str = ternary( p < 0.05, sprintf('p = %0.5f', p), 'ns' );

    ms = [nanmean(s0), nanmean(s1)];
    errs = [plotlabeled.nansem(s0), plotlabeled.nansem(s1)];
    [~, x] = plots.bars( axs(i), ms', PL{i, 2}, {p_str}, strrep(PL{i, 1}, '_', ' ') );
    plots.barerrs( axs(i), x, ms', errs' );
    set( axs(i), 'ylim', [-30e-3, 30e-3] );
  end
end

%%

t_sub = [ timeLabels >= 50 & timeLabels <= 350 ];
t_mean = squeeze( nanmean(coh(:, :, t_sub), 3) );
base_mask = rowmask( t_mean );

assert_ispair( t_mean, coh_labels );

plt_labels = coh_labels';
plt_mask = pipe( base_mask ...
  , @(m) find(plt_labels, 'acc_bla', m) ...
  , @(m) find(plt_labels, {'choice'}, m) ... 
  , @(m) find(plt_labels, {'cell-cluster1', 'cell-cluster2'}, m) ...
);

t_mean = t_mean(plt_mask, :);
plt_labels = plt_labels(plt_mask);

if ( 0 )
  sub_each = {'days', 'direction', 'lfp-channel', 'spk-channel', 'trialtype', 'cell-cluster'};
  [t_mean, plt_labels] = dsp3.summary_binary_op(...
    t_mean, plt_labels', sub_each, 'self-bottle', 'other-bottle', @minus, @nanmean );  
end

freq_win = [2, 60];
f_ind = find( f >= freq_win(1) & f <= freq_win(2) );

figure(1); clf;
pl = plotlabeled.make_common();

% pl.x = f(f_ind);
% axs = pl.lines( t_mean(:, f_ind), plt_labels, 'cell-cluster', {'trialtype', 'direction', 'outcome'} );

[I, id, C] = rowsets( 2, plt_labels, {'direction', 'cell-cluster'}, {'outcome'} );
[PI, PL] = plots.nest2( id, I, plots.cellstr_join(C) );
axs = plots.simple2( PI, PL, @(ax, inds, j, label) stdshade_cb(t_mean(:, f_ind), f(f_ind), inds, j, label) );
for i = 1:numel(PI)
  p0 = arrayfun( @(x) select(2, @() ttest2(t_mean(PI{i}{1}, x), t_mean(PI{i}{2}, x))), f_ind );
  sig = p0 < 0.05;
  scatter( axs(i), f(f_ind(sig)), max(get(axs(i), 'ylim')), 'k*' );
end
xlabel( axs, 'Hz' );
% shared_utils.plot.add_horizontal_lines( axs, 0 );
shared_utils.plot.set_ylims( axs, [-0.08, 0.08] );


%%

[I, id, C] = rowsets( 2, plt_labels, {'direction', 'cell-cluster'}, {'outcome'}, 'mask', plt_mask );
[PI, PL] = plots.nest2( id, I, plots.cellstr_join(C) );
axs = plots.simple2( PI, PL, @(ax, inds, j, label) stdshade_cb(f_mean, timeLabels, inds, j, label) );
shared_utils.plot.add_horizontal_lines( axs, 0 );
shared_utils.plot.set_xlims( axs, [-150, 600] );


%%

function p = ttest2p(a, b)
[~, p] = ttest2( a, b );
end

function h = stdshade_cb(data, timeLabels, inds, j, label)

colors = { 'r', 'b' };
dat = data(inds{j}, :);
h = stdshade( dat,.5,colors{j}, timeLabels, 3 );
% h = stdshade( dat,.5,plots.color(numel(inds), j), timeLabels, 3 );
set( h, 'displayname', label );

end

function [coh, coh_labels, f, t] = load_orig_spit(data_root)

sfcoh_p = fullfile( data_root, 'sfcoh-split-by-cell-type' );
k1_sfcoh = load( fullfile(sfcoh_p, 'all-trialtype_K1_SFC_Matrix.mat') );
k2_sfcoh = load( fullfile(sfcoh_p, 'all-trialtype_K2_SFC_Matrix.mat') );

%

[coh1, labels1] = linearize_sfcoh( k1_sfcoh.SFCMAT );
addsetcat( labels1, 'cell-cluster', 'cell-cluster1' );

[coh2, labels2] = linearize_sfcoh( k2_sfcoh.SFCMAT );
addsetcat( labels2, 'cell-cluster', 'cell-cluster2' );

coh = [ coh1; coh2 ];
coh_labels = append( labels1, labels2 );
assert_ispair( coh, coh_labels );

f = k1_sfcoh.loadStruct.coh_f;
t = cellfun(@mean,k1_sfcoh.loadStruct.binned_t)+k1_sfcoh.loadStruct.t(1);
assert( numel(f) == size(coh, 2) && numel(t) == size(coh, 3) );

%
t_ind = t >= -150 & t <= 650;

coh = coh(:, :, t_ind);
t = t(t_ind);

end

function [coh, coh_labels, f, t] = load_new_split(data_root, second_subset)

path2CohMat = fullfile( data_root, 'COH_MAT/COH_MATRIX_STDNORM.mat' );
coh = load( path2CohMat );

meanTimeBins = cellfun(@mean,coh.loadStruct.binned_t)+coh.loadStruct.t(1);
f = coh.loadStruct.coh_f;
t = meanTimeBins;
[coh, coh_labels, src_to_dst] = linearize_sfcoh( coh.averageMat );

t_ind = t >= -150 & t <= 650;
coh = coh(:, :, t_ind);
t = t(t_ind);

%
cluster_info = load( fullfile(data_root, 'cell-clusters/cell-clusters-reward-mag-121422.mat') );

s_info = shared_utils.io.fload( fullfile(data_root, 'COH_MAT/site_info.mat') );
copy_fields = { 'lfpChans', 'lfpDay', 'lfpRegs', 'spikeIdx', 'spkChan', 'spkReg' };
dst_fields = {'lfp-channel', 'days', 'lfp-region', 'unit_index', 'spk-channel', 'spk-region'};

s_info(:, 3) = cellfun( @(x) sprintf('lfp-%s', x), s_info(:, 3), 'un', 0 );
s_info(:, 6) = cellfun( @(x) sprintf('spk-%s', x), s_info(:, 6), 'un', 0 );
s_info = fcat.from( s_info, dst_fields );
s_info = s_info(src_to_dst);
join( coh_labels, s_info );

if ( second_subset )
  subset = load( fullfile(data_root, 'cell-clusters/backup.mat') );
  cluster_info = subset.cluster_info;
end

spike_info_I = bfw.find_combinations( coh_labels, table2array(cluster_info.spike_cluster_info)' );
for i = 1:numel(spike_info_I)
  addsetcat( coh_labels, 'cell-cluster', cluster_info.cell_cluster_labels{i}, spike_info_I{i} );
end

% keep_ind = find( coh_labels, {'acc_bla', 'spk-acc'} );
% coh_labels = coh_labels(keep_ind);
% coh = coh(keep_ind, :, :);

end

function sampled = find_cell_cluster_subsets_matching_n(coh_labels, base_mask, replace, up_sample)

unit_spec = {'spk-region', 'spk-channel', 'unit_index', 'days'};
cluster_I = findall( coh_labels, 'cell-cluster', findnone(coh_labels, '<cell-cluster>', base_mask) );
units_by_cluster = eachcell( @(x) combs(coh_labels, unit_spec, x), cluster_I );
ns = cellfun( @(x) size(x, 2), units_by_cluster );
min_n = min( ns );
max_n = max( ns );

if ( up_sample )
  assert( replace, 'Must allow replacement when upsampling' );
  sample_n = max_n;
else
  sample_n = min_n;
end

sampled = eachcell( @(x) x(:, randsample(size(x, 2), sample_n, replace)), units_by_cluster );
sampled = sort( cate1( eachcell(...
  @(x, m) cate1(bfw.find_combinations(coh_labels, x, m)), sampled, cluster_I)) );

end