function [pdc, pdc_ls, pdc_f] = run_pdc_pair(data, labels, target_pair, bi, varargin)

defaults = struct();
defaults.verbose = false;
defaults.n_freqs = 128;
defaults.is_orig_pdc_method = false;
defaults.asymp_pdc_params = {};
params = shared_utils.general.parsestruct( defaults, varargin );

%{

m(2, 1): 1 -> 2
m(1, 2): 2 -> 1

%}

nf = params.n_freqs;
nt = numel( target_pair{1} ); % num trials
assert( nt == numel( target_pair{2} ), 'Mismatching number of trials in pair.' );

if ( params.is_orig_pdc_method )
  % [nTrials, nSamples, nChans] = size(LFPData);
  chan0 = data(target_pair{1}, :);
  chan1 = data(target_pair{2}, :);
  X = cat( 3, chan0, chan1 );
  
  pdc = [];
  pdc_ls = fcat();
  pdc_f = [];
  
  for ti = 1:nt
    [tmp_pdc, ~, pdc_f] = partialDirectedCOH( X(ti, :, :) );
    pdc_mean = cellfun( @(x) nanmean(abs(tmp_pdc(:, :, :, x)), 4), bi, 'un', 0 );
    pdc_mean = cat( 4, pdc_mean{:} );
    dir0 = squeeze( pdc_mean(1, 2, :, :) );
    dir1 = squeeze( pdc_mean(2, 1, :, :) );
    
    ci = cellfun( @(x) x(ti), target_pair );
    append( pdc_ls, pair_labels(labels, ci) );
    
    if ( ti == 1 )
      pdc = nan( 2 * nt, numel(pdc_f), numel(bi) );
    end
    
    ti0 = (ti - 1) * 2;
    pdc(ti0 + 1, :, :) = dir0;
    pdc(ti0 + 2, :, :) = dir1;
  end
else
  pdc = cell( nt, 1 );
  pdc_ls = cell( nt, 1 );
  pdc_ft = cell( nt, 1 );
  pdc_f = [];

  for ti = 1:nt
    if ( params.verbose )
      fprintf( '\n trial %d of %d', ti, nt );
    end

    ci = cellfun( @(x) x(ti), target_pair );

    tmp_pdc = nan( 2, nf, numel(bi) );
    f = [];

    for bin = 1:numel(bi)
      X = data(ci, bi{bin});
      c = run_asymp_pdc( X, 'nFreqs', nf, params.asymp_pdc_params{:} );

      tmp_pdc(1, :, bin) = squeeze( c.pdc2(1, 2, :) );
      tmp_pdc(2, :, bin) = squeeze( c.pdc2(2, 1, :) );

      f = c.f;
    end

    pdc_ls{ti} = pair_labels( labels, ci );
    pdc{ti} = tmp_pdc;
    pdc_ft{ti} = f;
  end

  pdc = vertcat( pdc{:} );
  pdc_ls = vertcat( fcat, pdc_ls{:} );

  if ( ~isempty(pdc_ft) )
    pdc_f = pdc_ft{1};
  end
end

end

function res = pair_labels(labels, ci)

l_ci = labels(ci);
new_ls_dir0 = fcat.from( join_unique_labels_per_column(l_ci([1, 2])), getcats(l_ci) );
new_ls_dir1 = fcat.from( join_unique_labels_per_column(l_ci([2, 1])), getcats(l_ci) );
res = [ new_ls_dir0'; new_ls_dir1 ];

end