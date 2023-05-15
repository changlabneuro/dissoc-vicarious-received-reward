function [coh, coh_t, coh_f, coh_lfp, coh_spk, coh_uid, coh_idx, coh_bin] = SFC_wrapFunc(...
                                                                    lfp_data ...
                                                                  , lfp_labels ...
                                                                  , lfp_t ...
                                                                  , lfp_binned_t ...
                                                                  , spike_ts ...
                                                                  , spike_labels ...
                                                                  , event_times ...
                                                                  , event_labels ...
                                                                  , chronux_params)

%  lfp_cfg()                                                             
                                                              
assert_ispair( lfp_data, lfp_labels );
assert_ispair( spike_ts, spike_labels );
assert_ispair( event_times, event_labels );

% For each day, select the lfp data and units from that day. For each unit,
% compute a psth of its activity aligned to the events for that day. For
% each LFP-site, compute spike field coherence between the aligned spiking
% activity and the lfp data from the site.
[lfp_day_I, lfp_day_C] = findall( lfp_labels, 'days' );


coh = {};
phs = {};
coh_t = {};
coh_f = [];
coh_lfp = {};
coh_spk = {};
coh_uid =[];
coh_idx = {};
coh_bin = [];
cc = 0;% combo counter

for i = 1:numel(lfp_day_I)
  lfp_day_ind = lfp_day_I{i};
  % Find units from this day.
  unit_day_ind = find( spike_labels, lfp_day_C(:, i) );
  % Find events from this day.
  spike_event_ind = find( event_labels, lfp_day_C(:, i) );
  % Select the events from this day.
  spike_events = event_times(spike_event_ind);
  
  % Find all lfp sites on this day.
  [~, lfp_site_I] = ...
    keepeach( lfp_labels', {'channels', 'regions', 'sites'}, lfp_day_ind );
  
  % Loop over all combinations of sites and units.
  index_combinations = dsp3.numel_combvec( lfp_site_I, unit_day_ind );
  
  for j = 1:size(index_combinations, 2)
    fprintf( '\n %d of %d\t\t ', j, size(index_combinations, 2) );
    
    index_combination = index_combinations(:, j);
    lfp_site_ind = lfp_site_I{index_combination(1)};
    unit_ind = unit_day_ind(index_combination(2));
    
    assert( numel(spike_event_ind) == numel(lfp_site_ind) );
    spike_times = spike_ts{unit_ind};
    
    spk_info = spike_labels{unit_ind};
    spk_ID =  cellstr([spk_info(1, 5),  spk_info(1, 3),  spk_info(1, 12),spk_info(1, 17)]);
    
    lfp_info = lfp_labels{lfp_site_ind,:};
    
    lfp_ID = cellstr([lfp_info(:, 4),  lfp_info(:, 3),  lfp_info(:, 7),lfp_info(:, 9),lfp_info(:, 15), lfp_info(:, 1), lfp_info(:,11), ]); 
    
     % Increment the spike/lfp combo counter by 1
     cc = cc+1;
    
    % For each time window
    for k = 1:numel(lfp_binned_t)
      %fprintf( '\n\t %d of %d', k, numel(lfp_binned_t) );
      
      windowed_lfp = lfp_data(lfp_site_ind, lfp_binned_t{k});
      windowed_lfp_t = lfp_t(lfp_binned_t{k});
      windowed_spike_t_series = windowed_lfp_t * 1e-3; % to ms.
      
      t_min = spike_events + windowed_spike_t_series(1);
      t_max = spike_events + windowed_spike_t_series(end);
      
      is_spike_in_window = arrayfun( ...
        @(x, y) spike_times >= x & spike_times < y, t_min, t_max, 'un', 0 );
      % Select the spikes in this time window, and re-map the spike times
      % to [0, window_size]
      spikes_in_window = arrayfun( ...
          @(selector, t_min) columnize(spike_times(selector{1}) - t_min) ...
        , is_spike_in_window, t_min, 'un', 0 );
      
      chronux_spikes = struct( 'spikes', spikes_in_window(:) );
      
      % Call the chronux routine.
      [C, phi, ~, ~, ~, f] = coherencycpt( windowed_lfp', chronux_spikes', chronux_params );
      
     
      %fprintf( '\n Computed coherence; %d non-NaN observations.', sum(~isnan(C(:))) );
      %fprintf('%d.', sum(~isnan(C(:))) );
      fprintf('.');
      coh_f = f;
     
      coh{end+1, 1} = C;
      phs{end+1, 1} = phi;
      coh_lfp{end+1, 1} = lfp_ID;
      coh_spk{end+1, 1} = spk_ID;
      coh_bin(end+1, 1) = k;
      coh_idx{end+1, 1} = lfp_site_ind;
      coh_uid(end+1, 1) = cc;
      coh_t{end+1, 1} = windowed_lfp_t(1);
      
    end
    fprintf('\n');
  end
end

end