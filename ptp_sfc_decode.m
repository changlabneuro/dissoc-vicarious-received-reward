%% Load SFC data

path2CohMat = '/Volumes/external3/data/changlab/ptp-vicarious-reward/COH_MAT/COH_MATRIX_STDNORM.mat';

load(path2CohMat);

minPlotFreq = 10;
maxPlotFreq = 60;
minPlotMs = -250;
maxPlotMs = 750;

ab_f_ind = freqLabels > 10 & freqLabels < 20;
g_f_ind = freqLabels > 35 & freqLabels < 51;
t_ind = timeLabels > 50 & timeLabels < 400;

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

%%  decoding

ab_data = squeeze( nanmean(ptp_coh(:, ab_f_ind, t_ind), [2, 3]) );
g_data = squeeze( nanmean(ptp_coh(:, g_f_ind, t_ind), [2, 3]) );
src_labels = repset( addcat(ptp_labels', 'bands'), 'bands', {'alpha-beta', 'gamma'} );
src_data = [ ab_data; g_data ];
assert_ispair( src_data, src_labels );
src_mask = find( ~isnan(src_data) );

classify_each = { 'direction', 'trialtype', 'bands' };
[decode_labs, cls_I] = retaineach( src_labels, classify_each, src_mask );

nt = 1e2;

real_accs = cell( size(cls_I) );
null_accs = cell( size(cls_I) );
sig_ps = zeros( size(cls_I) );
parfor i = 1:numel(cls_I)
  fprintf( '\n %d of %d', i, numel(cls_I) );
  grp = categorical( src_labels, 'outcome', cls_I{i} );
  real_accs{i} = decode( src_data(cls_I{i}), grp, nt, false, true );  
  null_accs{i} = decode( src_data(cls_I{i}), grp, nt, true, true );
  sig_ps(i) = ranksum( real_accs{i}, null_accs{i} );
%   sig_ps(i) = sum( mean(real_accs{i}) < null_accs{i} ) / nt;
end

%%

ci_real = cellfun( @ci, real_accs, 'un', 0 );
ci_null = cellfun( @ci, null_accs, 'un', 0 );
sig_tf = cellfun( @(real, null) min(real) > max(null), ci_real, ci_null );

rep_labs = cate1( ...
  [{fcat}, arrayfun(@(i) repmat(decode_labs(i), nt), 1:numel(real_accs), 'un', 0)] );
dec_labs = repset( addcat(rep_labs, 'data-type'), 'data-type', {'real', 'null'} );
dec_acc = [ vertcat(real_accs{:}); vertcat(null_accs{:}) ];

p_mask = find( decode_labs, 'acc_bla' );
sig_ps(p_mask)
decode_labs(p_mask)

%%  box plots real vs null

figure( 1 );
clf;

plt_mask = find( dec_labs, {'acc_bla'} );
[I, id, C] = rowsets( 2, dec_labs, {'direction', 'bands', 'trialtype'}, 'data-type' ...
  , 'mask', plt_mask );
[PI, PL] = plots.nest2( id, I, plots.strip_underscore(plots.cellstr_join(C)) );

axs = plots.panels( numel(PI) );
for i = 1:size(PI, 1)
  [h, x] = plots.box( axs(i), dec_acc*1e2, PI{i}, PL{i, 2}, PL{i, 1} );
end

shared_utils.plot.match_ylims( axs );
hold( axs, 'on' );
shared_utils.plot.add_horizontal_lines( axs, 50 );
ylabel( axs(1), 'Accuracy (%)' )

%%

function accs = decode(src_data, grp, nt, shuff, n_match)

assert( numel(src_data) == numel(grp) );

accs = zeros( nt, 1 );
for i = 1:nt
  if ( shuff )
    grp = do_shuffle( grp );
  end
  
  if ( n_match )
    [~, ~, ic] = unique( grp );
    ic = groupi( ic );
    min_n = min( cellfun(@numel, ic) );
    ic = cate1( cellfun(@(x) keep_n(do_shuffle(x), min_n), ic, 'un', 0) );
    use_dat = src_data(ic);
    use_grp = grp(ic);
  else
    use_dat = src_data;
    use_grp = grp;
  end
  
  warn = warning( 'off', 'all' );
  part = cvpartition( use_grp, 'holdout', 0.3 );
  warning( warn );
  assert( numel(unique(use_grp(part.training))) == numel(unique(use_grp(part.test))) );
  mdl = fitcdiscr( use_dat(part.training), use_grp(part.training) );
  yp = predict( mdl, use_dat(part.test) );
  accs(i) = sum( yp == use_grp(part.test) ) / sum( part.test );
end

end

function x = do_shuffle(x)
x = x(randperm(numel(x)));
end

function x = keep_n(x, n)
x = x(1:n);
end

function c = ci(x)
sem = std(x)/sqrt(length(x));              
ts = tinv([0.025  0.975],length(x)-1);     
c = mean(x) + ts*sem;                      
end