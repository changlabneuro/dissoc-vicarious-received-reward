function nfft = get_nfft(N, pad)

nfft=max(2^(nextpow2(N)+pad),N); % number of points in fft of prolates

end