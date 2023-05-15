function [C, f] = sfcoherence_pair(lfp_channel, spike_ts, lfp_wi, spk_t_series_s, spk_bin_w_s, sua_event_ts, chronux_p)

assert( numel(spk_t_series_s) == numel(lfp_wi) );
assert( size(lfp_channel, 1) == numel(sua_event_ts) );
assert( isvector(spike_ts) && isa(spike_ts, 'double') );
assert( isscalar(spk_bin_w_s) );
assert( iscell(lfp_wi) && isa(spk_t_series_s, 'double') );

fpad = 0;
if ( isfield(chronux_p, 'pad') )
  fpad = chronux_p.pad;
end

f = get_fgrid( chronux_p.Fs, floor(spk_bin_w_s * 1e3), fpad );
C = nan( size(lfp_channel, 1), numel(f), numel(lfp_wi) );

for i = 1:numel(lfp_wi)
  spk_subset = bin_spikes( spike_ts ...
    , sua_event_ts + spk_t_series_s(i) - spk_bin_w_s * 0.5 ...
    , sua_event_ts + spk_t_series_s(i) + spk_bin_w_s * 0.5);
  lfp_subset = lfp_channel(:, lfp_wi{i});

  empty_spike_ts = arrayfun( @(x) isempty(x.times), spk_subset ) | isnan( sua_event_ts );
  if ( ~all(empty_spike_ts) )
    eval_data_a = lfp_subset(~empty_spike_ts, :);
    eval_data_b = spk_subset(~empty_spike_ts);
    
    [c, ~, ~, ~, ~, f] = coherencycpt( eval_data_a', eval_data_b', chronux_p );
    C(~empty_spike_ts, :, i) = c';
  end
end

end

function s = bin_spikes(spike_ts, event_start_ts, event_stop_ts)

assert( numel(event_start_ts) == numel(event_stop_ts) );

spikes = cell( numel(event_start_ts), 1 );
for i = 1:numel(event_start_ts)
  ib = spike_ts >= event_start_ts(i) & spike_ts < event_stop_ts(i);
  spikes{i} = struct( 'times', {columnize(sort(spike_ts(ib) - event_start_ts(i)))} );
end

s = vertcat( spikes{:} );

end