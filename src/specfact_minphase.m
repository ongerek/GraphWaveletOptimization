function s = specfact_minphase(tf, K)
%SPECFACT_MINPHASE Minimum-phase spectral factor of a nonnegative
%   symmetric Laurent polynomial via the cepstral (Kolmogorov) method.
%
%   s = SPECFACT_MINPHASE(tf, K):  tf = centered coefficients
%   [t_d ... t_1 t_0 t_1 ... t_d] of T(z) = T(1/z), T(e^{jw}) >= 0.
%   Returns s (length d+1) with  S(z) S(1/z) = T(z),  S minimum phase.
%   K = FFT size (default 2^16). T should be strictly positive for
%   full accuracy; tiny values are floored at 1e-14 * max(T).
tf = tf(:);
d  = (numel(tf) - 1) / 2;
if nargin < 2, K = 2^16; end
w  = 2*pi*(0:K-1).'/K;
Tw = real(exp(-1i * w * (-d:d)) * tf);
Tw = max(Tw, 1e-14 * max(Tw));
c  = real(ifft(log(Tw)));              % real cepstrum of T
slog          = zeros(K, 1);           % causal cepstrum of S
slog(1)       = c(1)/2;
slog(2:K/2)   = c(2:K/2);
slog(K/2+1)   = c(K/2+1)/2;
Sm = exp(fft(slog));                   % S_min on the grid
s  = real(ifft(Sm));
s  = s(1:d+1);
end
