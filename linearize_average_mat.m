function [ptp_coh, ptp_labels] = linearize_average_mat(averageMat, site_info, cons)

[ptp_coh, ptp_labels, src_to_dst] = linearize_sfcoh( averageMat );
site_ind = arrayfun( @(x) sprintf('site-%d', x), src_to_dst, 'un', 0 );

if ( ~isempty(site_info) )
  site_labs = site_info_to_labels( site_info );
  site_labs = rmcat( site_labs(src_to_dst) ...
    , {'lfp-region', 'spk-region', 'lfp-channel', 'spk-channel'} );
  addsetcat( site_labs, 'site_index', site_ind );

  join( ptp_labels, site_labs );
end

event_labs = fcat.from( cons.events.labels );

match_each = { 'days', 'monkeys' };
[match_I, match_C] = findall( event_labs, match_each );
for i = 1:numel(match_I)
  dst_ind = find( ptp_labels, match_C(1, i) );
  for j = 2:numel(match_each)
    addsetcat( ptp_labels, match_each{j}, match_C{j, i}, dst_ind );
  end
end

prune( ptp_labels );

end