function new_c = join_unique_labels_per_column(l)

c = cellstr( l );
new_c = cell( 1, size(c, 2) );
for i = 1:size(c, 2)
  new_c{i} = strjoin( unique(c(:, i), 'stable'), '_' );
end

end