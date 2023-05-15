function [all_accs, dec_labs] = ptp_decode_by_animal(ptp_coh, ptp_labels, freqLabels, timeLabels, mask)

assert_ispair( ptp_coh, ptp_labels );
if ( nargin < 5 )
  mask = rowmask( ptp_coh );
end

gfi = freqLabels > 35 & freqLabels < 51;
afi = freqLabels > 10 & freqLabels < 20;
ti = timeLabels > 50 & timeLabels < 400;
coh_gamma = squeeze( nanmean(ptp_coh(mask, gfi, ti, :), [2, 3]) );
coh_ab = squeeze( nanmean(ptp_coh(mask, afi, ti, :), [2, 3]) );
coh_win = [ coh_gamma; coh_ab ];
coh_labs = repset( addcat(ptp_labels(mask), 'band'), 'band', {'gamma', 'alpha_beta'} );
assert_ispair( coh_win, coh_labs );

dec_each = {'direction', 'trialtype', 'band', 'monkeys'};
[dec_I, dec_C] = findall( coh_labs, dec_each );
monks = combs( coh_labs, 'monkeys' );

its = 1e2;
ho = 0.3;

dec_labs = fcat();
all_accs = [];
for i = 1:numel(dec_I)
  fprintf( '\n %d of %d', i, numel(dec_I) );
  
  curr_monk = dec_C(end, i);
  self_ind = dec_I{i};
  other_ind = find( coh_labs, [dec_C(1:end-1, i)', setdiff(monks, curr_monk)] );
  monk_str = strjoin( [curr_monk, setdiff(monks, curr_monk)], '->' );
  
  grp_self = categorical( coh_labs, 'outcome', self_ind );
  grp_other = categorical( coh_labs, 'outcome', other_ind );
  
  accs = nan( its, 1 );
  for j = 1:its
    warn = warning( 'off', 'all' );
    part_self = cvpartition( grp_self, 'holdout', ho );
    part_other = cvpartition( grp_other, 'holdout', ho );
    warning( warn );
    
    si = self_ind(part_self.training);
    oi = other_ind(part_other.training);
    svm_res = fitcsvm( ...
      coh_win(si), grp_self(part_self.training), 'KernelFunction', 'rbf' );
    pv = predict( svm_res, coh_win(oi) );
    accs(j) = sum( pv == grp_other(part_other.training) ) / sum( part_other.training );
  end
  all_accs = [ all_accs; accs ];
  
  i0 = size( dec_labs, 1 ) + 1;
  append1( dec_labs, coh_labs, self_ind, numel(accs) );
  setcat( dec_labs, 'monkeys', monk_str, i0:size(dec_labs, 1) );
end

end