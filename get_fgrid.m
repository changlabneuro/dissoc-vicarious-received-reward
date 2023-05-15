function f = get_fgrid(Fs, N, pad)

if ( nargin < 3 )
  pad = -1;
end

f = getfgrid( Fs, get_nfft(N, pad), [0, Fs * 0.5] );

end