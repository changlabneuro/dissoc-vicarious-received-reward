function batchComputeSpikeFieldCohNewPad

%'pad'

dr = '/Volumes/external3/data/changlab/ptp-vicarious-reward';

lfp_file_p = '/Users/putnampt/Dropbox (ChangLab)/ptp/VicariousRewardLFP/reward-aligned';
outDir = '/Users/putnampt/Dropbox (ChangLab)/ptp/VicariousRewardLFP/paddedCOH';

lfp_file_p = fullfile( dr, 'reward_aligned_lfp/reward-aligned' );
outDir = fullfile( dr, 'remade-sfcoh-bigger-time-window' );

% lfp_cfg();
tc = 3047;
tc = 0;
lfp_dir = dir(fullfile(lfp_file_p, '*.mat')); %load( 'lfp_day__06152016_reward.mat' );

t_win = 300;  % ms
t_step = 50;  % ms

for f = 4:size(lfp_dir,1)
    
%     keep f lfp_dir tc
    
    
    lfp_file = load(fullfile(lfp_dir(f).folder, lfp_dir(f).name));
    sua_file = load( fullfile(dr, 'sua/dictator_game_SUAdata_pre.mat') );
    cons = load( 'trial_data.mat' );
    
    spike_info = dsp3_linearize_cc_sua_data( sua_file );
    [data, labels, categories, t] = signal_container_lfp_to_unformatted( lfp_file.lfp );
    
    lfp_dflts = dsp3.get_common_lfp_defaults();
    chronux_params = lfp_dflts.chronux_params;
    
    lfp_labels = fcat.from( labels, categories );
    spike_labels = spike_info.spike_labels';
    
    events = cons.consolidated.events;
    
    event_labels = fcat.from( events.labels );
    event_key = cons.consolidated.event_key;
    event_times = events.data(:, event_key('rwdOn'));
    spike_times = spike_info.spike_times;
    
    numDataSamples = size(data, 2);
    
    binned_t = shared_utils.vector.slidebin( 1:numDataSamples, t_win, t_step, true );
    
    numTimeBins = size(binned_t,2);
    
    [coh, coh_t, coh_f, coh_lfp, coh_spk, coh_uid, coh_idx, coh_bin] = SFC_wrapFunc( data, lfp_labels, t, binned_t ...
        , spike_times, spike_labels, event_times, event_labels, chronux_params );
    
    uniqCombos = unique(coh_uid);
    
    numFreq = size(coh_f,2);
    
    for i = 1:size(uniqCombos,1)
        
        cc = uniqCombos(i);
        ccIdxs = find(coh_uid == cc);
        
        if size(ccIdxs,1) ~=  numTimeBins, error('Number of indexs for spike/field does not match number of timepoints!');end
        
        ccSpikeInfo     = coh_spk{ccIdxs(1)};
        
        spkDay = ccSpikeInfo{1};
        spkChan = ccSpikeInfo{2};
        spkReg = ccSpikeInfo{3};
        spikeIdx = ccSpikeInfo{4};
        spkStr = sprintf('%s_%s_%s_%s',spkDay, spkChan, spkReg, spikeIdx);
        
        ccLFPInfo       = coh_lfp{ccIdxs(1)};
        
        numObservations = size(ccLFPInfo,1);
        
        lfpDays = ccLFPInfo(:,1);
        lfpChans = ccLFPInfo(:,2);
        lfpRegs = ccLFPInfo(:,7);
        
        if numel(unique(lfpDays))~=1 || numel(unique(lfpChans))~=1 || numel(unique(lfpRegs))~=1
            error('Incosistant Rows for LFP table')
            lfpStr = [];
        else
            lfpDay = lfpDays{1};
            lfpChans = lfpChans{1};
            lfpRegs = lfpRegs{1};
            lfpStr = sprintf('%s_%s_%s',lfpDay, lfpChans, lfpRegs);
        end
        
        cohMatrix = nan(numFreq, numTimeBins, numObservations);
        
        
        for ti = 1:numTimeBins
            
            cohMatrix(:,ti,:) = coh{ccIdxs(ti)};
        end
        
        tc = tc+1;
        
        saveStr =  sprintf('%d_%s_%s_%s_%s-%s-%s_COH.mat',tc, spkDay, spkChan, spkReg, spikeIdx, lfpChans, lfpRegs);
        savePath = fullfile(outDir,saveStr);
        
        save(savePath, 'cohMatrix', 'ccLFPInfo', 'ccSpikeInfo', 'spkDay', 'spkChan', 'spkReg', 'spikeIdx', 'lfpDay', 'lfpChans', 'lfpRegs', 'coh_f', 't', 'binned_t',  '-v7.3');
        
        
        
    end
end

% Inputs:
% Size of 'DATA'            =  4108 x 2150 ((num_observations x num_lfp_chans) x timePoints)
% Size of 'lfp_labels'      =  4108 x 15  ((num_observations x num_lfp_chans) x labels)
% Size of 't'               =     1 x 2150
% Size of 'binned_t'        =     1 x 41 --> (-1000:1:1149)
% Size of 'spike_times'     =   343 x 1  --> numRecordings x 1
% Size of 'spike_labels'    =   343 x 17 --> numRecordings x numRecordingAttributes
% Size of 'event_times'     = 61324 x 1
% Size of 'event_labels'    = 61324 x 13

% Outputs:
% Size of 'coh'     --> 246 x 1 cell, where each cell is size(data,1) x size(coh_f,2)
% Size of 'coh_t'   --> 246 x 1 cell (6 repeating series of -1K --> 1K steps of 50
% Size of 'coh_f'   --> 1 x 129 double, 0 --> 500 (hz?)

% 6 is the number of combinations of sites and units
% 41 time bins
% 246 = 41 time bins per each unique combo of site/unit (








