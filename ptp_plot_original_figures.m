clearvars
clc
%% Paths to COH matrices
% path2CohMat = '/Users/putnampt/Dropbox (ChangLab)/ptp/VicariousRewardLFP/COH_MAT/COH_MATRIX_STDNORM.mat';
% path2CohMat = '/Volumes/external3/ptp-vicarious-reward/COH_MAT/COH_MATRIX_STDNORM.mat';
path2CohMat = '/Volumes/external3/data/changlab/ptp-vicarious-reward/COH_MAT/COH_MATRIX_STDNORM.mat';

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

ttest_bin_by_bin = false;

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

acc_bla = true;

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

shared_utils.plot.set_clims( axH([1, 2, 5, 6]), [-.014, 0.013] );

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