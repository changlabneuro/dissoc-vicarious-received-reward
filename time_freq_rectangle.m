function h = time_freq_rectangle(ax, t, f, ts, fs)

assert( numel(ts) == 2 && numel(fs) == 2 );

xt = get( ax, 'xtick' );
yt = get( ax, 'ytick' );
assert( numel(xt) == numel(t) && numel(yt) == numel(f) );

x0 = interp1( t, xt, ts(1) );
x1 = interp1( t, xt, ts(2) );

y1 = min( max(f), interp1(flip(f), yt, fs(1)) );
y0 = interp1( flip(f), yt, fs(2) );

h = rectangle( ax, 'position', [x0, y0, x1-x0, y1-y0] );

end