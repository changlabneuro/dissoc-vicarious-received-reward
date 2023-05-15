function [cons_coh, cons_labels, sign_coh, sign_labels] = ptp_sign_consistency(ptp_coh, ptp_labels, vic_gt_rcv, mask)

assert_ispair( ptp_coh, ptp_labels );

mu_each = {'direction', 'monkeys', 'outcome', 'trialtype'};
[mean_labels, I] = retaineach( ptp_labels, mu_each, mask );
mus = bfw.row_nanmean( ptp_coh, I );

[sign_labels, I] = retaineach( mean_labels, setdiff(mu_each, 'outcome') );
sign_coh = nan( [numel(I), size(ptp_coh, 2:3)] );
for i = 1:numel(I)
  s_ind = find( mean_labels, 'self-bottle', I{i} );
  o_ind = find( mean_labels, 'other-bottle', I{i} );
  
  if ( vic_gt_rcv )
    sign_coh(i, :, :) = mus(o_ind, :, :) > mus(s_ind, :, :);
  else
    sign_coh(i, :, :) = mus(o_ind, :, :) < mus(s_ind, :, :);
  end
end

[cons_labels, I] = retaineach( sign_labels, setdiff(mu_each, 'monkeys') );
cons_coh = nan( [numel(I), size(ptp_coh, 2:3)] );
for i = 1:numel(I)
  h_ind = find( sign_labels, 'hitch', I{i} );
  k_ind = find( sign_labels, 'kuro', I{i} );
  assert( isscalar(h_ind) && isscalar(k_ind) );
  
  h_coh = sign_coh(h_ind, :, :);
  k_coh = sign_coh(k_ind, :, :);
  
  match = h_coh == k_coh;
  match_pos = match & h_coh == 1;
  match_neg = match & h_coh == 0;
  
  cons_coh(i, match_pos) = 1;
  cons_coh(i, match_neg) = 2;
  cons_coh(i, h_coh == 1 & ~match) = 3;
  cons_coh(i, k_coh == 1 & ~match) = 4;
end

end