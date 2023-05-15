function [tot_coh, tot_coh_labels] = do_norm_contrast(coh, coh_labels, I, across, really_contrast)

if ( nargin < 5 )
  really_contrast = true;
end

assert_ispair( coh, coh_labels );

if ( isempty(coh) )
  tot_coh = [];
  tot_coh_labels = fcat();
  return
end

copy_labels_sn = coh_labels';
copy_labels_on = coh_labels';
setcat( copy_labels_sn, 'outcomes', 'self-none' );
setcat( copy_labels_on, 'outcomes', 'other-none' );

tot_coh = cell( size(I) );
tot_coh_labels = cell( size(I) );

for i = 1:numel(I)
  fprintf( '\n%d of %d', i, numel(I) );
  
  si = I{i};
  tt_I = findall( coh_labels, {'trialtypes', 'administration'}, si );
  
  tmp_coh = [];
  tmp_labels = fcat();
  
  for j = 1:numel(tt_I)
    ti = tt_I{j};
    s_ind = find( coh_labels, {'self'}, ti );
    o_ind = find( coh_labels, {'other'}, ti );
    n_ind = find( coh_labels, {'none'}, ti );
    
    coh_s = coh(s_ind, :, :);
    coh_n = coh(n_ind, :, :);
    coh_o = coh(o_ind, :, :);
    
    coh_sn = nanmean( [coh_s; coh_n], 1 );
    coh_on = nanmean( [coh_o; coh_n], 1 );
    
    rest_each = setdiff( getcats(coh_labels), csunion(across, 'outcomes') );
    rest_I = findall( coh_labels, rest_each, ti );
    for k = 1:numel(rest_I)
      ri = rest_I{k};
      
      s_ind = find( coh_labels, 'self', ri );
      o_ind = find( coh_labels, 'other', ri );
      n_ind = find( coh_labels, 'none', ri );
      
      if ( ~really_contrast )  % no contrast, just subtract mean
        coh_s = (nanmean(coh(s_ind, :, :), 1) - coh_sn);
        coh_o = (nanmean(coh(o_ind, :, :), 1) - coh_on);
        coh_n = (nanmean(coh(n_ind, :, :), 1) - coh_on);

        append1( tmp_labels, coh_labels, s_ind );
        append1( tmp_labels, coh_labels, o_ind );
        append1( tmp_labels, coh_labels, n_ind );
        
        if ( isempty(n_ind) )
          coh_n = [];
        end
        if ( isempty(s_ind) )
          coh_s = [];
        end
        if ( isempty(o_ind) )
          coh_o = [];
        end
        
        tmp_coh = [ tmp_coh; coh_s; coh_o; coh_n ];
      else
        if ( numel(s_ind) + numel(n_ind) == 1 )
          coh_s = nan( 1, size(coh, 2), size(coh, 3) );
        elseif ( numel(s_ind) + numel(n_ind) == 0 )
          coh_s = [];
        else
          coh_s = (nanmean(coh(s_ind, :, :), 1) - coh_sn) - ...
            (nanmean(coh(n_ind, :, :), 1) - coh_sn);
        end
        if ( numel(o_ind) + numel(n_ind) == 1 )
          coh_o = nan( 1, size(coh, 2), size(coh, 3) );
        elseif ( numel(o_ind) + numel(n_ind) == 0 )
          coh_o = [];
        else
          coh_o = (nanmean(coh(o_ind, :, :), 1) - coh_on) - ...
            (nanmean(coh(n_ind, :, :), 1) - coh_on);
        end

        append1( tmp_labels, copy_labels_sn, [s_ind; n_ind] );      
        append1( tmp_labels, copy_labels_on, [o_ind; n_ind] );
        tmp_coh = [ tmp_coh; coh_s; coh_o ];
      end
      
      assert_ispair( tmp_coh, tmp_labels );
    end
  end
  
  tot_coh{i} = tmp_coh;
  tot_coh_labels{i} = tmp_labels;
end

tot_coh = vertcat( tot_coh{:} );
tot_coh_labels = vertcat( fcat, tot_coh_labels{:} );

end