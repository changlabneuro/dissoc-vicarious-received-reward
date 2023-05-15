function ind = n_match(labels, sep_I, within, across)

for idx = 1:numel(sep_I)
  I = findall( labels, within, sep_I{idx} );
  
  ns = cellfun( @numel, I, 'un', 0 );
  assert( isequal(ns{:}), 'Expect matching numbers of elements within subsets.' );
  
  cts = cellfun( @(x) cellfun(@numel, findall(labels, across, x)), I, 'un', 0 );
  assert( isequal(cts{:}), 'Expect equal frequencies of conditions within subset.' );
  
  for i = 1:numel(I)
    xi = findall( labels, across, I{i} );
    if ( i == 1 )
      % sample a random set of trials
      min_n = min( cellfun(@numel, xi) );
      si = cellfun( @(x) sort(randsample(numel(x), min_n)), xi, 'un', 0 );
    end
    % use `si` to sample the same set of trials for each subset `I`
    I{i} = cate1( cellfun(@(x, i) x(i), xi, si, 'un', 0) );
  end
  
  sep_I{idx} = vertcat( I{:} );
end

ind = sort( vertcat(sep_I{:}) );

end