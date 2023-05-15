function [pair_sets, pairs] = find_pairs(labels, of, specific_to, mask)

%   FIND_PAIRS -- Find pairs of 

if ( nargin < 4 )
  mask = rowmask( labels );
end

region_subsets = findeach( labels, of, mask );
region_pairi = bfw.pair_combination_indices( numel(region_subsets) );

pair_sets = {};
pairs = {};
for i = 1:size(region_pairi, 1)
  [i0, c0] = findeach( labels, specific_to, region_subsets{region_pairi(i, 1)} );
  [i1, c1] = findeach( labels, specific_to, region_subsets{region_pairi(i, 2)} );
  for j = 1:numel(i0)
    for k = 1:numel(i1)
      pair_sets(end+1, :) = [ i0(j), i1(k) ];
      pairs(end+1, :) = [ c0(j), c1(k) ];
    end
  end
end

end