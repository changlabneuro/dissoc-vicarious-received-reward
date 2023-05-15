clear all
clc
 
binSizesMs = [150, 50];
timeWindowMs = 1000;
normMethod = 'STDNORM';
numKlustas = 2;
klustaStr = 'everything';

data_root = '/Volumes/external3/data/changlab/ptp-vicarious-reward';
 
% saveDir = fullfile(pwd, 'plots');
% addpath(genpath(fullfile(pwd,'categorical-master')))
% addpath(genpath(fullfile(pwd,'global-master')))
% addpath(genpath(fullfile(pwd,'dsp-master')))
% addpath(genpath(fullfile(pwd,'dsp3-master')))
% addpath(genpath(fullfile(pwd,'shared_utils-master')))
 
SFCTablePath = '/Users/putnampt/dev/pastdoc/vicariousRewardLFP'; %'H:\dev\postdoc\vicariousRewardLFP'
 
 
% desiredTimeWinS = [.05 .250];
desiredTimeWinS = [.05 .350];
 
sua_file    = load( fullfile(data_root, 'sua/dictator_game_SUAdata_pre.mat') );
cons      = load( fullfile(data_root, 'trial_data/trial_data.mat') );
spike_cluster_info   = dsp3_linearize_cc_sua_data( sua_file );
spike_labels  = spike_cluster_info.spike_labels';
events     = cons.consolidated.events;
event_labels  = fcat.from( events.labels );
event_key    = cons.consolidated.event_key;
event_times   = events.data(:, event_key('rwdOn'));
spike_times   = spike_cluster_info.spike_times;

assert_ispair( spike_times, spike_labels );
assert_ispair( event_times, event_labels );
 
nUnits = size(spike_times,1);
 
[spike_day_I, spike_day_C] = findall( spike_labels, 'days' );
[event_day_I, event_day_C] = findall( event_labels, 'days' );
 
ExpVarMat = nan(nUnits, 4);
 
unitUUIDs = cell(nUnits,1);
 
codingMetrics = struct();
 
pcaMat = nan(nUnits, 0);
 
 
 
%  n = numel(uniqueResponseValues);
%   
%   r = zeros(size(uniqueResponseValues));
%   
%   for j = 1:n
%    thisUniqueResponseIdxs = find(not(cellfun('isempty', strfind(withinLimitsResponsesValues, uniqueResponseValues{j}))));
%     
%    r(j) = mean(withinLimitsPredictorValues(thisUniqueResponseIdxs));
% 
%   end
%   
%   dos = ( n - ( sum(r) / max(r) ) ) / (n-1);

%%
 
for u = 1:nUnits
  fprintf( '\n %d of %d', u, nUnits );
   
  %%
  
  spike_day    = spike_labels(u, 'days');
  spike_reg    = spike_labels(u, 'regions');
  spike_chan   = spike_labels(u, 'channels');
  spike_unit   = spike_labels(u, 'unit_index');
   
  
   
  unitUUIDs{u} = sprintf('%s_%s_%s_%s', spike_day{1}, spike_chan{1},spike_reg{1}, spike_unit{1});
   
   codingMetrics(u).UUID = unitUUIDs{u};
    codingMetrics(u).Day = spike_day;
    codingMetrics(u).SpikeRegion = spike_reg;
    codingMetrics(u).SpikeChannel = spike_chan;
    codingMetrics(u).Unit = spike_unit;
     
    
   
  event_idxs   = find(event_labels, spike_day);
   
  spike_event_labels = event_labels(event_idxs);
   
  unit_event_times  = event_times(event_idxs);
  unit_spike_times  = spike_times{u};
   
  unit_event_administration  = event_labels(event_idxs, 'administration');
  unit_event_blocks      = event_labels(event_idxs, 'blocks');
  unit_event_contexts     = event_labels(event_idxs, 'contexts');
  unit_event_days       = event_labels(event_idxs, 'days');
  unit_event_drugs      = event_labels(event_idxs, 'drugs');
  unit_event_magnitudes    = event_labels(event_idxs, 'magnitudes');
  unit_event_monkeys     = event_labels(event_idxs, 'monkeys');
  unit_event_outcomes     = event_labels(event_idxs, 'outcomes');
  unit_event_recipients    = event_labels(event_idxs, 'recipients');
  unit_event_sessionID    = event_labels(event_idxs, 'session_ids');
  unit_event_sessions     = event_labels(event_idxs, 'sessions');
  unit_event_trialnumber   = event_labels(event_idxs, 'trials');
  unit_event_trialtype    = event_labels(event_idxs, 'trialtypes');
   
  preIdxs   = find(strcmpi(unit_event_administration, 'pre'));
  selfIdxs  = find(strcmpi(unit_event_outcomes, 'self'));
  otherdxs  = find(strcmpi(unit_event_outcomes, 'other'));
  bothIdxs  = find(strcmpi(unit_event_outcomes, 'both'));
  noneIdxs  = find(strcmpi(unit_event_outcomes, 'none'));
  cuedIdxs  = find(strcmpi(unit_event_trialtype, 'cued'));
  choiceIdxs = find(strcmpi(unit_event_trialtype, 'choice'));
   
  highIdxs  = find(strcmpi(unit_event_magnitudes, 'high'));
  mediumIdxs = find(strcmpi(unit_event_magnitudes, 'medium'));
  lowIdxs   = find(strcmpi(unit_event_magnitudes, 'low'));
   
   
  cuedIdxs = intersect(cuedIdxs, preIdxs);
  choiceIdxs = intersect(choiceIdxs, preIdxs);
   
   
   
  type_conditionIdxs{1} = intersect(selfIdxs, choiceIdxs);
  type_conditionIdxs{2} = intersect(otherdxs, choiceIdxs);
  type_conditionIdxs{3} = intersect(bothIdxs, choiceIdxs);
  type_conditionIdxs{4} = intersect(noneIdxs, choiceIdxs);
  type_conditionIdxs{5} = intersect(selfIdxs, cuedIdxs);
  type_conditionIdxs{6} = intersect(otherdxs, cuedIdxs);
  type_conditionIdxs{7} = intersect(bothIdxs, cuedIdxs);
  type_conditionIdxs{8} = intersect(noneIdxs, cuedIdxs);
   
  type_conditionLabels{1} = 'Self (Choice)';
  type_conditionLabels{2} = 'Other (Choice)';
  type_conditionLabels{3} = 'Both (Choice)';
  type_conditionLabels{4} = 'None (Choice)';
  type_conditionLabels{5} = 'Self (Forced)';
  type_conditionLabels{6} = 'Other (Forced)';
  type_conditionLabels{7} = 'Both (Forced)';
  type_conditionLabels{8} = 'None (Forced)';
   
   
  type_conditionIdxs{9} = intersect(preIdxs,selfIdxs);
  type_conditionIdxs{10} = intersect(preIdxs,otherdxs);
  type_conditionIdxs{11} = intersect(preIdxs,bothIdxs);
  type_conditionIdxs{12} = intersect(preIdxs,noneIdxs);
   
   
  type_conditionLabels{9} = 'Self (All)';
  type_conditionLabels{10} = 'Other (All)';
  type_conditionLabels{11} = 'Both (All)';
  type_conditionLabels{12} = 'None (All)';
   
   
  type_conditionIdxs{13} = intersect(preIdxs,highIdxs);
  type_conditionIdxs{14} = intersect(preIdxs,cuedIdxs);
  type_conditionIdxs{15} = intersect(preIdxs,lowIdxs);
  type_conditionLabels{13} = 'High reward (All)';
  type_conditionLabels{14} = 'Medium Reward (All)';
  type_conditionLabels{15} = 'Low Reward (All)';
   
   
   
   
  type_conditionIdxs{16} = intersect(preIdxs,choiceIdxs);
  type_conditionIdxs{17} = intersect(preIdxs,cuedIdxs);
  type_conditionLabels{16} = 'Choice (All)';
  type_conditionLabels{17} = 'Cued (All)';
  
%   for h = 1:numel(type_conditionIdxs)
%     [type_evokedFRs, binMeanS] = getEvokedFR(unit_spike_times,unit_event_times(type_conditionIdxs{h}), binSizesMs(:), timeWindowMs);  
%   end
   
  [type_evokedFRs, binMeanS] = getEvokedFRs(unit_spike_times,unit_event_times, type_conditionIdxs, timeWindowMs, binSizesMs);
  [type_conditionFRs, ~] = getConditionalFRs(unit_spike_times,unit_event_times, type_conditionIdxs, timeWindowMs, binSizesMs);
   
   norm_conditionFRs = (type_conditionFRs - min(min(type_conditionFRs))) / ( max(max(type_conditionFRs)) - min(min(type_conditionFRs)) );
   
  desiredTimeIdxs = find( binMeanS > desiredTimeWinS(1) & binMeanS < desiredTimeWinS(2));
   
 
  pcaMat(u, 1:size(type_conditionFRs,1)) = mean(norm_conditionFRs(:,desiredTimeIdxs),2);
   
  desiredTimeIdxs = find( binMeanS > desiredTimeWinS(1) & binMeanS < desiredTimeWinS(2));
  %% Outcome / Choice Trials
  y = [];
  y = [mean(type_evokedFRs{1}(desiredTimeIdxs,:),1)'; ...
    mean(type_evokedFRs{2}(desiredTimeIdxs,:),1)'; ...
    mean(type_evokedFRs{3}(desiredTimeIdxs,:),1)'; ...
    mean(type_evokedFRs{4}(desiredTimeIdxs,:),1)'];
   
  grp = {};
  grp = vertcat(repmat({'SELF'}, size(type_evokedFRs{1},2),1),...
    repmat({'OTHR'}, size(type_evokedFRs{2},2),1));
  grp = vertcat(grp, repmat({'BoTH'}, size(type_evokedFRs{3},2),1));
  grp = vertcat(grp,repmat({'NONE'}, size(type_evokedFRs{4},2),1));
   
   
  [p,tbl,stats,terms] = anovan(y,{grp}, 'display','off');
   
  EV_ChoiceOutcome = (tbl{2,2}/tbl{4,2});
  codingMetrics(u).EV_ChoiceOutcome = EV_ChoiceOutcome;
   
  ExpVarMat(u,1) = EV_ChoiceOutcome;
   
  %% Outcome / Cued Trials
  y = [];
  y = [mean(type_evokedFRs{5}(desiredTimeIdxs,:),1)'; ...
    mean(type_evokedFRs{6}(desiredTimeIdxs,:),1)'; ...
    mean(type_evokedFRs{7}(desiredTimeIdxs,:),1)'; ...
    mean(type_evokedFRs{8}(desiredTimeIdxs,:),1)'];
   
  grp = {};
  grp = vertcat(repmat({'SELF'}, size(type_evokedFRs{5},2),1),...
    repmat({'OTHR'}, size(type_evokedFRs{6},2),1));
  grp = vertcat(grp, repmat({'BoTH'}, size(type_evokedFRs{7},2),1));
  grp = vertcat(grp,repmat({'NONE'}, size(type_evokedFRs{8},2),1));
   
   
  [p,tbl,stats,terms] = anovan(y,{grp}, 'display','off');
   
  EV_CuedOutcome = (tbl{2,2}/tbl{4,2});
   codingMetrics(u).EV_CuedOutcome = EV_CuedOutcome;
  ExpVarMat(u,2) = EV_CuedOutcome;
  %% Reward Magnitude / All Trials
    y = [];
  y = [mean(type_evokedFRs{13}(desiredTimeIdxs,:),1)'; ...
    mean(type_evokedFRs{14}(desiredTimeIdxs,:),1)'; ...
    mean(type_evokedFRs{15}(desiredTimeIdxs,:),1)'];
   
  grp = {};
  grp = vertcat(repmat({'HI'}, size(type_evokedFRs{13},2),1),...
    repmat({'ME'}, size(type_evokedFRs{14},2),1));
  grp = vertcat(grp, repmat({'LO'}, size(type_evokedFRs{15},2),1));
 
   
   
  [p,tbl,stats,terms] = anovan(y,{grp}, 'display','off');
   
  EV_RewardMag = (tbl{2,2}/tbl{4,2});
   codingMetrics(u).EV_RewardMag = EV_RewardMag;
  ExpVarMat(u,3) = EV_RewardMag;
   
  %% Cued Vs Choice Trials
   y = [];
  y = [mean(type_evokedFRs{16}(desiredTimeIdxs,:),1)'; ...
    mean(type_evokedFRs{17}(desiredTimeIdxs,:),1)'];
   
  grp = {};
  grp = vertcat(repmat({'Choice'}, size(type_evokedFRs{16},2),1),...
    repmat({'Forced'}, size(type_evokedFRs{17},2),1));
 
   
   
  [p,tbl,stats,terms] = anovan(y,{grp}, 'display','off');
   
  EV_TrialType = (tbl{2,2}/tbl{4,2});
    codingMetrics(u).EV_TrialType = EV_TrialType;
  ExpVarMat(u,4) = EV_TrialType;
   
  %% Self vs Other choice trials
   y = [];
  y = [mean(type_evokedFRs{1}(desiredTimeIdxs,:),1)'; ...
    mean(type_evokedFRs{2}(desiredTimeIdxs,:),1)'];
   
  grp = {};
  grp = vertcat(repmat({'Choice'}, size(type_evokedFRs{1},2),1),...
    repmat({'Forced'}, size(type_evokedFRs{2},2),1));
 
   
   
  [p,tbl,stats,terms] = anovan(y,{grp}, 'display','off');
   
  EV_SelfVsOtherChoice = (tbl{2,2}/tbl{4,2});
   codingMetrics(u).EV_SelfVsOtherChoice = EV_SelfVsOtherChoice;
  ExpVarMat(u,5) = EV_TrialType;
   
  %% Choice-Outcome DOS
   
  r = [mean(mean(type_evokedFRs{1}(desiredTimeIdxs,:),1)); ...
    mean(mean(type_evokedFRs{2}(desiredTimeIdxs,:),1)); ...
    mean(mean(type_evokedFRs{3}(desiredTimeIdxs,:),1)); ...
    mean(mean(type_evokedFRs{4}(desiredTimeIdxs,:),1))];
   
  n = size(r,1);
   
  dos = ( n - ( sum(r) / max(r) ) ) / (n-1);
   
  codingMetrics(u).DOS_ChoiceOutcome = dos;
   
  %%
end
 
[idx,C] = kmeans(pcaMat,numKlustas);
[coeff,score,latent,tsquared,explained,mu] = pca(pcaMat);

%%

metrics = SUACoding{:, 6:11};
[idx, C] = kmeans( metrics, 2 );
%%  remake clusters

SUACoding = struct2table(codingMetrics);

target_cells = strcmp( SUACoding.SpikeRegion, 'acc' );
split_by = { 'EV_TrialType', 'EV_RewardMag', 'DOS_ChoiceOutcome'  };
metrics = table2array( SUACoding(:, split_by) );
metrics = metrics(target_cells, :);
[~, score] = pca( metrics );
[idx, C] = kmeans( score, 2 );

spike_cluster_info = SUACoding(target_cells, {'SpikeChannel', 'Unit', 'Day'});
cell_cluster_labels = repmat( {'cell-cluster1'}, numel(idx), 1 );
cell_cluster_labels(idx == 2) = { 'cell-cluster2' };

clf; figure(1);
hold on;
scatter( score(idx == 1, 1), score(idx == 1, 2), 'r*' );
scatter( score(idx == 2, 1), score(idx == 2, 2), 'b*' );

if ( 0 )
  save( fullfile(data_root, 'cell-clusters/cell-clusters-reward-mag-121422.mat'), 'idx', 'spike_cluster_info', 'cell_cluster_labels', 'metrics', 'score' );
end

%%

cluster_info = load( fullfile(data_root, 'cell-clusters/cell-clusters-reward-mag-121422.mat') );
score = cluster_info.score;
idx = cluster_info.idx;

%%

clf();
histogram( score(idx == 1, 1), 10 ); 
hold on;
histogram( score(idx == 2, 1), 10 );
xlim( [-0.55, 0.55] );
figure( 1 );
 
%%
 
% codingMetrics(:).PCA1Score = score(:,1);
% codingMetrics.PCA2Score = score(:,2);
% codingMetrics.PCA3Score = score(:,3);
 
%%
SUACoding = struct2table(codingMetrics);
 
SUACoding.PCA1Score = score(:,1);
SUACoding.PCA2Score = score(:,2);
SUACoding.PCA3Score = score(:,3);
save('SUACoding3.mat', 'SUACoding');
 
% %%
% [idx,C] = kmeans(pcaMat,numKlustas);
% [coeff,score,latent,tsquared,explained,mu] = pca(pcaMat);
%
% fig1H = figure(1);
% clf
% subplot(8,4,1:20)
% hold on
% gscatter(score(:,1),score(:,2),idx,'bgm')
% xlabel('PC1');
% ylabel('PC2');
% legend('Cluster 1','Cluster 2')
%
% subplot(8,4,21:22)
% hold on
% histogram( score((idx==1),1))
% histogram( score((idx==2),1))
% xlabel('PC1');
%
% subplot(8,4,23:24)
% hold on
% histogram( score((idx==1),2))
% histogram( score((idx==2),2))
% xlabel('PC2');
%
%
% subplot(8,4,25:26)
% hold on
% histogram( score((idx==1),3))
% histogram( score((idx==2),3))
% xlabel('PC3');
% %
% % subplot(8,4,27:28)
% % hold on
% % histogram( score((idx==1),4))
% % histogram( score((idx==2),4))
% % xlabel('PC4');
%
% %
% % subplot(8,4,29:30)
% % hold on
% % histogram( score((idx==1),5))
% % histogram( score((idx==2),5))
% % xlabel('PC5');
% %
% %
% % subplot(8,4,31:32)
% % hold on
% % histogram( score((idx==1),6))
% % histogram( score((idx==2),6))
% % xlabel('PC6');
%
% sgtitle(klustaStr);
% makeLandscapePDF(fig1H, sprintf('%s_Clusters.pdf', klustaStr));
%
%
%
%
%
%
%
%