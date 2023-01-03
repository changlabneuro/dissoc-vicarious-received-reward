ff_mats = shared_utils.io.findmat( '/Volumes/external3/data/changlab/ptp-vicarious-reward/f2fgrang' );

% ff_mats = ff_mats(1);

for i = 1:numel(ff_mats)
  shared_utils.general.progress( i, numel(ff_mats) );
  
  ff_granger = shared_utils.io.fload( ff_mats{i} );
  
  for j = 1:numel(ff_granger.Results)    
    chan = sprintf( '%s_%s', char(ff_granger.blaChan{j}), char(ff_granger.accChan{j}) );
    addsetcat( ff_granger.Results(j).labels, 'bla_acc_channels', chan );
  end
  
  labels = vertcat( ff_granger.Results.labels );
  granger = vertcat( ff_granger.Results.data );
  assert( isequal(unique(cellfun(@numel, findall(labels, {'region', 'bla_acc_channels', 'outcomes'}))), 1) );
  
  mask = pipe( rowmask(labels) ...
    , @(m) findnone(labels, 'errors', m) ...
  );

  granger = granger(mask, :, :);
  labels = keep( labels, mask );
  
  assert_ispair( granger, labels );
  f = ff_granger.Results(1).f;
  t = ff_granger.Results(1).t;
  keep_f = f >= 10 & f <= 100;
  f = f(keep_f);
  granger = granger(:, keep_f, :);
  
  dst_p = fullfile( fileparts(fileparts(ff_mats{i})), 'f2fgrang-smaller' ...
    , shared_utils.io.filenames(ff_mats{i}, true) );
  save( dst_p, 'granger', 'labels', 'f', 't' );
end

