function [all_accs, dec_labs] = ptp_decode_per_animal(ptp_coh, ptp_labels, freqLabels, timeLabels, mask)

assert_ispair( ptp_coh, ptp_labels );

if ( nargin < 5 )
  mask = rowmask( ptp_coh );
end

gfi = freqLabels > 35 & freqLabels < 51;
afi = freqLabels > 10 & freqLabels < 20;
% ti = find( timeLabels > 50 & timeLabels < 400 );
ti = true( size(timeLabels) );

coh_gamma = squeeze( nanmean(ptp_coh(mask, gfi, :, :), 2) );
coh_ab = squeeze( nanmean(ptp_coh(mask, afi, :, :), 2) );
coh_win = [ coh_gamma; coh_ab ];
coh_labs = repset( addcat(ptp_labels(mask), 'band'), 'band', {'gamma', 'alpha_beta'} );
assert_ispair( coh_win, coh_labs );

dec_each = {'direction', 'trialtype', 'band', 'monkeys'};
[dec_I, dec_C] = findall( coh_labs, dec_each );

its = 1e2;
ho = 0.3;

dec_labs = fcat();
all_accs = [];

for i = 1:numel(dec_I)
  fprintf( '\n %d of %d', i, numel(dec_I) );
  
  curr_monk = dec_C(end, i);
  self_ind = dec_I{i};
  grp_self = categorical( coh_labs, 'outcome', self_ind );
  
  accs = nan( its, numel(ti) );
  for idx = 1:numel(ti)
    fprintf( '\n\t %d of %d', idx, numel(ti) );
    for j = 1:its      
      warn = warning( 'off', 'all' );
      part_self = cvpartition( grp_self, 'holdout', ho );
      warning( warn );

      si = self_ind(part_self.training);
      oi = self_ind(part_self.test);
      svm_res = fitcdiscr( coh_win(si, idx), grp_self(part_self.training) );
      pv = predict( svm_res, coh_win(oi, idx) );
      accs(j, idx) = sum( pv == grp_self(part_self.test) ) / sum( part_self.test );
    end
  end
  all_accs = [ all_accs; accs ];
  append1( dec_labs, coh_labs, self_ind, size(accs, 1) );
end

end