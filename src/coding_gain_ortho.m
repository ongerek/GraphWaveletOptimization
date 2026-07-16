function [GdB, info] = coding_gain_ortho(designer, r, L)
%CODING_GAIN_ORTHO L-level dyadic subband coding gain, orthonormal FB.
%
%   [GdB, info] = CODING_GAIN_ORTHO(designer, r, L)
%   designer : function handle  @(r) -> h0  called at every level with the
%              autocorrelation of the current lowpass signal (pass
%              @(r) h0_fixed to use one fixed filter at all levels);
%   r        : input autocorrelation (r(1) = r[0]);   L : levels.
%
%   G_L = sigma_x^2 / [ (sigma_{aL}^2)^{2^-L} * prod_j (sigma_{dj}^2)^{2^-j} ]
%
%   Uses the exact subband-variance recursion
%       sigma_a^2 = sum_k f[k] r[k],   r_a[m] = sum_k f[k] r[2m - k],
%   with f = conv(h0, fliplr(h0)) (h1 variance follows from the
%   orthonormal energy split sigma_d^2 = 2 r[0] - sigma_a^2).
r = r(:);
sigx2 = r(1);
sd2 = zeros(L,1);
for j = 1:L
    h0 = designer(r);
    n  = numel(h0);
    f  = conv(h0(:).', fliplr(h0(:).'));
    c0 = n;
    sa2    = f(c0)*r(1) + 2 * (f(c0+1:c0+n-1) * r(2:n));
    sd2(j) = 2*r(1) - sa2;
    Lr = floor((numel(r) - n) / 2);
    rn = zeros(Lr+1, 1);
    ks = -(n-1):(n-1);
    fv = f(c0 + ks);
    for m = 0:Lr
        rn(m+1) = fv * r(abs(2*m - ks) + 1);
    end
    r = rn;
end
saL2 = r(1);
den  = saL2^(2^-L) * prod(sd2 .^ (2.^-(1:L)'));
GdB  = 10*log10(sigx2 / den);
info = struct('sd2', sd2, 'saL2', saL2, 'sigx2', sigx2);
end
