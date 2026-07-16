function r = acf_from_psd(Sfun, K)
%ACF_FROM_PSD Autocorrelation from a PSD:  r[k] = (1/2pi) int S e^{jwk} dw.
%   r = ACF_FROM_PSD(Sfun, K): Sfun = handle S(w) on [0, 2pi), vectorized;
%   K = FFT size (default 2^15). Returns r(1)=r[0], r(2)=r[1], ...
if nargin < 2, K = 2^15; end
w = 2*pi*(0:K-1).'/K;
r = real(ifft(Sfun(w)));
r = r(1:K/2);
end
