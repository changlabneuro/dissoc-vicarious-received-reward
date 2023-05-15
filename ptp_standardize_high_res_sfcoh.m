dr = '/Volumes/external3/data/changlab/ptp-vicarious-reward';

is_high_res = false;

if ( is_high_res )
  mats = shared_utils.io.findmat( fullfile(dr, 'high-res-sfcoh') );
  % mats = mats(2);
  dst_p = fullfile( dr, 'high-res-sfcoh-contrast' );
else
  mats = shared_utils.io.findmat( fullfile(dr, 'remade-sfcoh') );
  dst_p = fullfile( dr, 'remade-sfcoh-contrast' );
end

%

for i = 1:numel(mats)
  fprintf( '\n %d of %d', i, numel(mats) );
  
  coh_file = load( mats{i} );
  
  %%
  
  t_mask = coh_file.t >= -0.05;
  coh = coh_file.coh(:, :, t_mask);
  labels = coh_file.labels';
  
  mask = pipe( rowmask(labels) ...
    , @(m) findnone(labels, 'errors', m) ...
    , @(m) find(labels, 'pre', m) ...
  );
  
  %%
  
  norm_each = { 'unit_id', 'lfp-channel', 'administration', 'trialtypes' };
  [contrast_coh, contrast_labs] = do_simple_norm_contrast( coh, labels, norm_each, mask );
  
  coh_file.coh = contrast_coh;
  coh_file.labels = contrast_labs;
  coh_file.t = coh_file.t(t_mask);
  
  if ( 1 )
    save_fname = shared_utils.io.filenames( mats{i}, true );
    save( fullfile(dst_p, save_fname), 'coh_file' );
  end  
end