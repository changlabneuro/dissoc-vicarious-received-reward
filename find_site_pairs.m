function m = find_site_pairs(site_pairs, labels)

m = [];

for i = 1:numel(site_pairs.days)
  if ( isempty(find(labels, site_pairs.days{i})) )
    continue;
  end
  
  target_pairs = site_pairs.channels{i};
  dir0 = fcat.strjoin( target_pairs', '_' );
  dir1 = fcat.strjoin( fliplr(target_pairs)', '_' );
  dirs = [ dir0'; dir1' ];
  
  for j = 1:size(dirs, 1)
    m = union( m, find(labels, [dirs(j, :), site_pairs.days(i)]) );
  end
end

end