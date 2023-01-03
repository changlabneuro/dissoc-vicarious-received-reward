conf = dsp3.set_dataroot( '/Volumes/external3/data/changlab/dsp3' );
outs = dsp3_find_iti_looks( 'config', conf, 'is_parallel', true, 'verbose', true );
cons = dsp3.get_consolidated_data( conf );

%%

addsetcat( outs.labels, 'session_ids', fullcat(fcat.from(cons.events.labels), 'session_ids') );

pl2_dir = 'H:\SIGNALS\raw';
[sesh_I, sesh_C] = findall( outs.labels, {'days', 'session_ids'} );

for i = 1:numel(sesh_I)
  match_ind = strcmp( cons.pl2_info.sessions, sesh_C{2, i} );
  assert( nnz(match_ind) == 1 );
  pl2_filename = cons.pl2_info.files{match_ind};
  
  pl2_file = shared_utils.io.find( pl2_dir, pl2_filename, true );
	assert( ~isempty(pl2_file) );
  
  chans_by_region = cons.pl2_info.channel_map(pl2_filename);
  tot_num_channels = sum( arrayfun(@(x) numel(x.channels), chans_by_region) );
  tot_lfp_data = [];
  lfp_index = 1;
  tot_lfp_labels = {};
  
  for j = 1:numel(chans_by_region)
    reg = chans_by_region(j);
    for k = 1:numel(reg.channels)
      chan_str = dsp3.channel_n_to_str( 'WB', reg.channels(k) );
      chan_data = PL2Ad( pl2_file, chan_str );
      
      if ( isempty(tot_lfp_data) )
        tot_lfp_data = nan( tot_num_channels, numel(chan_data.Values) );
      end
      
      tot_lfp_data(lfp_index, :) = chan_data.Values;
      tot_lfp_labels(end+1, :) = { reg.region, chan_str };
    end
  end
  
  aligned_lfp = dsp3.aligned_lfp_custom_events( ...
    files, picto_event_times, picto_event_labels, cons );
end

%%

tot_mag_size_data = {};
tot_mag_size_labels = {};

for i = 1:size(sub_fs,1)
  
fprintf( '\n%d of %d', i, size(sub_fs, 1) );

loadStruct = load(sub_fs{i});

if strcmp(loadStruct.spkReg, loadStruct.lfpRegs)
  continue;
end

lfp_ls = make_labels( loadStruct );

end

function lfp_ls = make_labels(loadStruct)

lfp_labs = loadStruct.ccLFPInfo;
lfp_labs(:, 2) = cellfun( @(x) sprintf('lfp-%s', x), lfp_labs(:, 2), 'un', 0 );
lfp_cols = { 'days', 'lfp-channels', 'magnitudes' ...
  , 'outcomes', 'trialtypes', 'administration', 'lfp-regions' };
lfp_ls = fcat.from( lfp_labs, lfp_cols );
assert( rows(lfp_ls) == size(loadStruct.cohMatrix, 3) );
spk_labs = loadStruct.ccSpikeInfo;
spk_labs(:, 3) = cellfun( @(x) sprintf('spk-%s', x), spk_labs(:, 2), 'un', 0 );
spk_ls = fcat.from( spk_labs, {'days', 'spk-channels', 'spk-regions', 'unit_index'} );
join( lfp_ls, spk_ls );

end