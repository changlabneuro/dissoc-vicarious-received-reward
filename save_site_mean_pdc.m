%%

dr = '/gpfs/milgram/project/chang/CHANG_LAB/naf3/Data/ptp-vicarious-reward';
src_p = fullfile( dr, 'pdc' );
dst_p = fullfile( dr, 'pdc-site-mean-n-match' );
shared_utils.io.require_dir( dst_p );

site_each = {'sites', 'channels', 'regions', 'outcomes', 'trialtypes', 'administration'};

do_n_match = true;
match_each = setdiff( site_each, {'outcomes'} );

src_mats = shared_utils.io.findmat( src_p );
for i = 1:numel(src_mats)
  dst_mat = fullfile( dst_p, shared_utils.io.filenames(src_mats{i}, true) );
  if ( exist(dst_mat, 'file') )
    continue;
  end
  
  mat = load( src_mats{i} );
  
  if ( do_n_match )
    I = findall( mat.pdc_ls, match_each );
    for j = 1:numel(I)
      outs = findall( mat.pdc_ls, 'outcomes', I{j} );
      min_n = min( cellfun(@numel, outs) );
      outs = cellfun( @(x) randsample(x, min_n), outs, 'un', 0 );
      I{j} = sort( vertcat(outs{:}) );
    end
    I = vertcat( I{:} );
    keep( mat.pdc_ls, I );
    mat.pdc = mat.pdc(I, :, :);
  end
  
  [mat.pdc_ls, I] = keepeach( mat.pdc_ls, site_each );
  mat.pdc = bfw.row_nanmean( real(mat.pdc), I );
  save( dst_mat, '-struct', 'mat' );
end


