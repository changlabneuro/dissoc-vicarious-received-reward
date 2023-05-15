dr = '/Volumes/external3/data/changlab/ptp-vicarious-reward';

lfp_mats = shared_utils.io.findmat( fullfile(dr, 'reward_aligned_lfp/reward-aligned') );
fnames = shared_utils.io.filenames( lfp_mats );
day_labs = cellfun( @(x) x(numel('lfp_')+1:numel('lfp_')+13), fnames, 'un', 0 );

site_info = site_info_to_labels( shared_utils.io.fload(fullfile(dr, 'COH_MAT/site_info.mat')) );

sua = load( fullfile(dr, 'sua/dictator_game_SUAdata_pre.mat') );
sua = dsp3_linearize_cc_sua_data( sua );

cons = shared_utils.io.fload( 'trial_data.mat' );

event_labs = fcat.from( cons.events.labels );
reward_ts = cons.events.data(:, cons.event_key('rwdOn'));

%%

is_high_res = false;
skip_existing = false;

if ( is_high_res )
  dst_p = fullfile( dr, 'high-res-sfcoh' );
  lfp_stp_size_ms = 10;
else
  dst_p = fullfile( dr, 'remade-sfcoh' );
  lfp_stp_size_ms = 50;
end
% lfp_win_size_ms = 150;
lfp_win_size_ms = 300;
lfp_t_series_ms = -1e3:lfp_stp_size_ms:1e3;

p = dsp3.make.defaults.coherence;
chronux_p = p.chronux_params;

if ( is_high_res )
  chronux_p.pad = 2;
end

min_t_ms = -50;
max_t_ms = 500;

subset_t = lfp_t_series_ms >= min_t_ms & lfp_t_series_ms <= max_t_ms;
spk_t_series_s = lfp_t_series_ms(subset_t) / 1e3;
spk_bin_w_s = lfp_win_size_ms / 1e3;

[day_I, day_C] = findall( sua.spike_labels, {'days', 'session_ids'} );

%

for i = 1:numel(day_I)
  dst_fname = sprintf( '%s.mat', strjoin(day_C(:, i), '_') );
  dst_fp = fullfile( dst_p, dst_fname );
  if ( skip_existing && exist(dst_fp, 'file') )
    fprintf( 'Skipping "%s"; already exists.\n', dst_fp );
    continue;    
  end
  
  %%
  
  fprintf( '\n %d of %d', i, numel(day_I) );
  
  di = strcmp( day_labs, day_C(1, i) );
  assert( sum(di) == 1 );
  lfp_file = shared_utils.io.fload( lfp_mats{di} );
  
  sua_event_ind = find( event_labs, day_C(:, i) );
  sua_event_ts = reward_ts(sua_event_ind);
  sua_event_ts(sua_event_ts == 0) = nan;
    
  lfp = lfp_file.data;
  lfp_labs = fcat.from( lfp_file.labels );
  
  if ( 1 )
    chn_cts = cellfun(@numel, findall(lfp_labs, 'channels'), 'un', 0);
    assert( isequal(chn_cts{:}) );
    if ( chn_cts{1} ~= numel(sua_event_ind) )
      fprintf( '\n Trial count mismatch.' );
      
      pre_ind = find( lfp_labs, 'pre' );
      chn_cts = cellfun(@numel, findall(lfp_labs, 'channels', pre_ind), 'un', 0);      
      if ( ~isempty(chn_cts) && isequal(chn_cts{:}, numel(sua_event_ind)) )
        fprintf( ' Using pre indices.\n' );
        pre_mask = pre_ind;
      else
        pre_ind = find( lfp_labs, 'session__1' );        
        chn_cts = cellfun(@numel, findall(lfp_labs, 'channels', pre_ind), 'un', 0);
        assert( ~isempty(chn_cts) && isequal(chn_cts{:}, numel(sua_event_ind)) );
        fprintf( ' Using session__1 indices.\n' );
        pre_mask = pre_ind;
      end
      
      lfp = lfp(pre_mask, :);
      lfp_labs = prune( lfp_labs(pre_mask) );
    else
      fprintf( '\n Trial count match.' );
    end
  end
  
  lfp_wi = shared_utils.vector.slidebin( ...
    1:size(lfp, 2), lfp_win_size_ms, lfp_stp_size_ms, true );
  lfp_wi = lfp_wi(subset_t);
  
  %%
  
  [~, src_pairs] = findeach( site_info, {'lfp-channel', 'unit_index', 'spk-region'} ...
    , find(site_info, day_C(1, i)) );
  
  try 
    %%
    [C, F, T, clfp_labs, cspk_labs] = run_pairs( sua, lfp, lfp_labs, src_pairs, day_C(:, i) ...
      , lfp_wi, spk_t_series_s, spk_bin_w_s, sua_event_ts, chronux_p );
  catch err
    rethrow( err );    
    warning( err.message );
    continue;
  end
  
  coh_labs = repmat( event_labs(sua_event_ind), size(src_pairs, 1) );  
  coh_labs = make_labels( coh_labs, clfp_labs', cspk_labs' );
  assert_ispair( C, coh_labs );
  coh_file = struct( 'coh', C, 'f', F, 't', T, 'labels', coh_labs );
  
  if ( 1 )
    fprintf( '\nSaving "%s".', dst_fp );
    save( dst_fp, '-v7.3', '-struct', 'coh_file' );
    fprintf( ' ... Done\n' );
  end
end

%%

function coh_labs = make_labels(event_labs, clfp_labs, cspk_labs)

lfp_regs = eachcell( @(x) sprintf('lfp-%s', x), clfp_labs(:, 'regions') );
spk_regs = eachcell( @(x) sprintf('spk-%s', x), cspk_labs(:, 'regions') );

setcat( clfp_labs, 'regions', lfp_regs );
setcat( cspk_labs, 'regions', spk_regs );

renamecat( clfp_labs, 'regions', 'lfp-region' );
renamecat( clfp_labs, 'channels', 'lfp-channel' );

renamecat( cspk_labs, 'regions', 'spk-region' );
renamecat( cspk_labs, 'channels', 'spk-channel' );

coh_labs = join( event_labs, clfp_labs, cspk_labs );
  
end

function [oc, of, ot, olfp_labs, ospk_labs] = run_pairs(...
  sua, lfp, lfp_labs, src_pairs, day ...
  , lfp_wi, spk_t_series_s, spk_bin_w_s, sua_event_ts, chronux_p)

max_f = 80;

C = cell( size(src_pairs, 1), 1 );
F = cell( size(C) );
T = cell( size(C) );

olfp_labs = cell( size(C) );
ospk_labs = cell( size(C) );

parfor j = 1:size(src_pairs, 1)
  fprintf( '\n\t %d of %d', j, size(src_pairs, 1) );

  lfp_ind = find( lfp_labs, src_pairs{j, 1} );
  spk_ind = find( sua.spike_labels, [day(:)', strrep(src_pairs(j, 2:end), 'spk-', '')] );

  assert( numel(lfp_ind) == numel(sua_event_ts) );
  assert( numel(spk_ind) == 1 );

  lfp_chan = lfp(lfp_ind, :);
  spk_ts = sua.spike_times{spk_ind};
  
%   try
%     assert( abs(min(spk_ts) - min(sua_event_ts(sua_event_ts > 0))) < 30 );
%   catch err
%     d = 10;
%   end

  [C{j}, f] = sfcoherence_pair( ...
    lfp_chan, spk_ts, lfp_wi, spk_t_series_s, spk_bin_w_s, sua_event_ts, chronux_p );
  keep_f = f <= max_f;
  f = f(keep_f);
  C{j} = C{j}(:, keep_f, :);
  F{j} = f;
  T{j} = spk_t_series_s;
  olfp_labs{j} = lfp_labs(lfp_ind);
  ospk_labs{j} = repmat( sua.spike_labels(spk_ind), numel(lfp_ind) );
end

oc = vertcat( C{:} );
of = F{1};
ot = T{1};
olfp_labs = vertcat( fcat, olfp_labs{:} );
ospk_labs = vertcat( fcat, ospk_labs{:} );

prune( olfp_labs );
prune( ospk_labs );

end