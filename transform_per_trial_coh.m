function [coh, labels, f, t] = transform_per_trial_coh(coh, labels, f, t)

mask = pipe( rowmask(labels) ...
  , @(m) findnone(labels, 'errors', m) ...
  , @(m) find(labels, 'pre', m) ...
);

site_each = getcats( labels );
[labels, coh_I] = keepeach( labels', site_each, mask );
coh = bfw.row_nanmean( coh, coh_I );

if ( ~isempty(labels) )
  addsetcat( labels, 'site_uuid', shared_utils.general.uuid );
end

end