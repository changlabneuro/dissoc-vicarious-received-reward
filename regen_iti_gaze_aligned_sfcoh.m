lfp_mats = shared_utils.io.findmat( '/Volumes/external3/data/changlab/ptp-vicarious-reward/iti_gaze_aligned_lfp' );

dst_p = fullfile( fileparts(fileparts(lfp_mats{1})), 'iti_gaze_aligned_coh' );
shared_utils.io.require_dir( dst_p );

% sua = load( '/Volumes/external3/data/changlab/ptp-vicarious-reward/sua/dictator_game_SUAdata_pre.mat' );
sua = load( 'dictator_game_SUAdata_pre.mat' );
sua = dsp3_linearize_cc_sua_data( sua );

renamecat( sua.spike_labels, 'channels', 'channel' );
renamecat( sua.spike_labels, 'regions', 'region' );
renamecat( sua.spike_labels, 'unit_id', 'unit_uuid' );

for i = 1:numel(lfp_mats)
  shared_utils.general.progress( i, numel(lfp_mats) );
  
  try
    lfp = load( lfp_mats{i} );

    non_finite = any( ~isfinite(lfp.tot_lfp_data), 2 );
    if ( all(non_finite) )
      fprintf( '\n Skipping' );
      continue;
    end

    lfp.tot_lfp_data(non_finite, :) = [];
    keep( lfp.tot_lfp_labels, find(~non_finite) );

    num_events = unique( cellfun('prodofsize', findall(lfp.tot_lfp_labels, 'channels')) );
    assert( numel(num_events) == 1 );

    if ( 1 )
      event_ts = lfp.tot_plex_events(1:num_events);
    else
      event_ts = rand( num_events, 1 );
    end

    lfp_file = struct();
    lfp_file.data = lfp.tot_lfp_data;
    lfp_file.labels = lfp.tot_lfp_labels;
    lfp_file.sample_rate = 1e3;
    lfp_file.params = struct( 'window_size', 0.15, 'step_size', 0.05, 'min_t', -0.5, 'max_t', 0.5 );
    lfp_file.src_filename = char( combs(lfp_file.labels, 'pl2_filename') );
    lfp_file.event_times = event_ts;
    renamecat( lfp.tot_lfp_labels, 'channels', 'channel' );

    out = dsp3.make.sfcoherence( struct('lfp', lfp_file), 'lfp', sua.spike_times, sua.spike_labels, 'verbose', true );
    save( fullfile(dst_p, shared_utils.io.filenames(lfp_mats{i}, true)), 'out', '-v7.3' );
  catch err
    warning( err.message );
  end
end

%%

