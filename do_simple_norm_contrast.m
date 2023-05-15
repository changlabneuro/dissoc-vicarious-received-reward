function [contrast_coh, contrast_labs] = do_simple_norm_contrast(coh, labels, norm_each, varargin)

assert_ispair( coh, labels );

[norm_labs, norm_I] = retaineach( labels, norm_each, varargin{:} );
  
contrast_labs = fcat();
contrast_coh = nan( numel(norm_I)*2, size(coh, 2), size(coh, 3) );

for j = 1:numel(norm_I)
  ni = norm_I{j};
  s_ind = find( labels, 'self', ni );
  o_ind = find( labels, 'other', ni );
  n_ind = find( labels, 'none', ni );

  if ( 1 )
    sn = nanmean( coh(s_ind, :, :), 1 ) - nanmean( coh(n_ind, :, :), 1 );
    on = nanmean( coh(o_ind, :, :), 1 ) - nanmean( coh(n_ind, :, :), 1 );
  else
    sn_b = nanmean( coh([s_ind; n_ind], :, :), 1 );
    sn = nanmean( coh(s_ind, :, :) - sn_b, 1 ) - nanmean( coh(n_ind, :, :) - sn_b, 1 );

    on_b = nanmean( coh([o_ind; n_ind], :, :), 1 );
    on = nanmean( coh(o_ind, :, :) - on_b, 1 ) - nanmean( coh(n_ind, :, :) - on_b, 1 );
  end

  sn_labs = retaineach( norm_labs, {}, j );
  setcat( sn_labs, 'outcomes', 'self-none' );

  on_labs = retaineach( norm_labs, {}, j );
  setcat( on_labs, 'outcomes', 'other-none' );
  append( contrast_labs, sn_labs );
  append( contrast_labs, on_labs );

  contrast_coh((j-1)*2+1, :, :) = sn;
  contrast_coh((j-1)*2+2, :, :) = on;
end

assert_ispair( contrast_coh, contrast_labs );

end