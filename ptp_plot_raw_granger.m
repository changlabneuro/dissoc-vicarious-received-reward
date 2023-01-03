%%
clear all
clc

data_root = '/Volumes/external3/data/changlab/ptp-vicarious-reward';

load( fullfile(data_root, 'granger_prepro/GrangerPrePro.mat') );
 
%%

ROI_MinFreq = 8;
ROI_MaxFreq = 20;

ROI_MinFreq = 35;
ROI_MaxFreq = 51;

%
 
% norm_self_acc_bla_mat = self_acc_bla_mat-none_acc_bla_mat;
% norm_self_bla_acc_mat = self_bla_acc_mat-none_bla_acc_mat;
%  
% norm_other_acc_bla_mat = other_acc_bla_mat-none_acc_bla_mat;
% norm_other_bla_acc_mat = other_bla_acc_mat-none_bla_acc_mat;
 
tsMin = 0;
tsMax = 0.1;
 
ROI_TimeIdxs = find(tLabs >= ROI_MinMs & tLabs <= ROI_MaxMs);
ROI_FreqIdxs = find(fLabs >= ROI_MinFreq & fLabs <= ROI_MaxFreq);
 
ROI_TimeRange = tLabs(ROI_TimeIdxs(end))-tLabs(ROI_TimeIdxs(1));
ROI_FreqRange = fLabs(ROI_FreqIdxs(end))-fLabs(ROI_FreqIdxs(1));
%%
figH = figure(1);
clf
 
axesH(1) = subplot(2,1,1)
title({'Self vs Other vs Bottle trials', 'ACC ---> BLA'});
[self_acc_bla_H,other_acc_bla_H, none_acc_bla_H,self_acc_bla_Vals, other_acc_bla_Vals,none_acc_bla_Vals] = plotTripartSpectroTimeSeriesROI(axesH(1) , self_acc_bla_mat,other_acc_bla_mat,none_acc_bla_mat, ROI_FreqIdxs, tLabs);
xlim([-500 500]);
xlabel('Time (ms)');
legend([self_acc_bla_H other_acc_bla_H, none_acc_bla_H],{'Self','Other', 'Bottle'},'Location', 'best')
ylim([tsMin tsMax]);
ylim([-0.0 .1])
vline(0, 'k-')
axesH(2) = subplot(2,1,2)
 
title({'Self vs Other vs Bottle trials', 'BLA ---> ACC'});
[self_bla_acc_H,other_bla_acc_H,none_bla_acc_H, self_bla_acc_Vals, other_bla_acc_Vals, none_bla_acc_Vals] = plotTripartSpectroTimeSeriesROI(axesH(2) , self_bla_acc_mat,other_bla_acc_mat,none_bla_acc_mat, ROI_FreqIdxs, tLabs);
xlim([-500 500]);
xlabel('Time (ms)');
legend([self_bla_acc_H other_bla_acc_H, none_bla_acc_H],{'Self','Other','None'},'Location', 'best')
ylim([tsMin tsMax]);
ylim([-0.0 .1])
vline(0, 'k-')

%
sval_sets = { self_acc_bla_Vals, self_bla_acc_Vals };
nval_sets = { none_acc_bla_Vals, none_bla_acc_Vals };
oval_sets = { other_acc_bla_Vals, other_bla_acc_Vals };

for i = 1:numel(sval_sets)

svals = sval_sets{i};
nvals = nval_sets{i};
ovals = oval_sets{i};
ax = axesH(i);

[~, ps_sn] = arrayfun( @(x) ttest2(svals(:, x), nvals(:, x)), 1:size(ovals, 2) );
[~, ps_on] = arrayfun( @(x) ttest2(ovals(:, x), nvals(:, x)), 1:size(ovals, 2) );
[~, ps_so] = arrayfun( @(x) ttest2(svals(:, x), ovals(:, x)), 1:size(ovals, 2) );
sig_sn = find( ps_sn < 0.05 );
sig_on = find( ps_on < 0.05 );
sig_so = find( ps_so < 0.05 );

lim_diff = diff( get(gca, 'ylim') );
h_sn = plot( ax, tLabs(sig_sn), repmat(max(get(gca, 'ylim')), numel(sig_sn), 1), 'k*' );
set( h_sn, 'color', 'r' );
h_on = plot( ax, tLabs(sig_on), repmat(max(get(gca, 'ylim')) - lim_diff * 0.1, numel(sig_on), 1), 'k*' );
set( h_on, 'color', 'b' );
h_so = plot( ax, tLabs(sig_so), repmat(max(get(gca, 'ylim')) - lim_diff * 0.2, numel(sig_so), 1), 'k*' );
set( h_so, 'color', 'k' );

end

%%

B = [8, 20]; 
G = [35, 51];

sab = self_acc_bla_mat;
sba = self_bla_acc_mat;
oab = other_acc_bla_mat;
oba = other_bla_acc_mat;
nab = none_acc_bla_mat;
nba = none_bla_acc_mat;

l = @(o, dir, x) [repmat({o, dir}, size(x, 1), 1), arrayfun(@(x) sprintf('site-%d', x), (1:size(x, 1))', 'un', 0)];
mu = @(g, f) squeeze(nanmean(real(g(:, fLabs >= f(1) & fLabs <= f(2), :)), 2));

granger_vals = [...
    mu(sab, B); mu(oab, B); mu(nab, B) ...
    ; mu(sba, B); mu(oba, B); mu(nba, B) ...
    ; mu(sab, G); mu(oab, G); mu(nab, G) ...
    ; mu(sba, G); mu(oba, G); mu(nba, G) ...
];
granger_labs = [...
    l('self', 'acc->bla', self_acc_bla_Vals); 
    l('other', 'acc->bla', other_acc_bla_Vals); 
    l('none', 'acc->bla', none_acc_bla_Vals);
    l('self', 'bla->acc', self_bla_acc_Vals);
    l('other', 'bla->acc', other_bla_acc_Vals); 
    l('none', 'bla->acc', none_bla_acc_Vals) ...
];
granger_labs = fcat.from( granger_labs, {'outcome', 'direction', 'site'} );
granger_labs = repset( addcat(granger_labs, 'band'), 'band', {'beta', 'gamma'} );

%%

mean_granger = nanmean( granger_vals(:, tLabs >= 50 & tLabs <= 350), 2 );

spec = setdiff( getcats(granger_labs), 'outcome' );
[sng, sn_labs] = dsp3.summary_binary_op( ...
  mean_granger, granger_labs', spec, 'self', 'none', @minus, @identity );
setcat( sn_labs, 'outcome', 'self-none' );
[ong, on_labs] = dsp3.summary_binary_op( ...
  mean_granger, granger_labs', spec, 'other', 'none', @minus, @identity );
setcat( on_labs, 'outcome', 'other-none' );
contrast_granger = [ sng; ong ];
contrast_labs = [ sn_labs'; on_labs ];

pl = plotlabeled.make_common();
[axs, inds] = pl.bar( contrast_granger, contrast_labs, 'outcome', 'direction', 'band' );

for i = 1:numel(inds)
  ind = inds{i};
  for j = 1:size(ind, 1)
    a = contrast_granger(ind{j, 1});
    b = contrast_granger(ind{j, 2});
    [~, p] = ttest2( a, b );
    if ( p < 0.05 )
      hold( axs(i), 'on' );
      plot( axs(i), j, max(get(axs(i), 'ylim')) - 4e-3, 'k*' );
    end
  end
end

%%

mean_granger = nanmean( granger_vals(:, tLabs >= 50 & tLabs <= 350), 2 );
pl = plotlabeled.make_common();
pl.x_order = { 'self', 'other', 'none' };
[axs, inds] = pl.bar( mean_granger, granger_labs, 'outcome', 'direction', 'band' );

pthreshs = [0.05, 0.01, 0.001, 0.0001];
stars = { '*', '**', '***', '****' };

for i = 1:numel(inds)
  ind = inds{i};
  for j = 1:size(ind, 1)
    a = granger_vals(ind{j, 1});
    b = granger_vals(ind{j, 2});
    [~, p] = ttest2( a, b );
    
    if ( p < pthreshs(end) )
      text( axs(i), j, max(get(axs(i), 'ylim')) - 4e-3, stars{end} );
    else
      for k = 2:numel(pthreshs)
        if ( p < pthreshs(k-1) && p >= pthreshs(k) )        
          text( axs(i), j, max(get(axs(i), 'ylim')) - 4e-3, stars{k} );
          break;
        end
      end
    end
  end
end

%%
makeLandscapePDF(figH, 'AlphaBeta_GrangerDirectionalityWithBottle.pdf')