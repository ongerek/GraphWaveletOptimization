function h0 = db_filter(p)
%DB_FILTER Minimum-phase Daubechies orthonormal filter with p vanishing
%   moments (length 2p), generated from the maxflat half-band product
%   filter via cepstral spectral factorization (no toolbox needed).
y  = [-1 2 -1]/4;                      % 1 - B(z)
tf = zeros(1, 2*(p-1)+1); yk = 1;
for k = 0:p-1
    coef = 2 * nchoosek(p-1+k, k);
    dd   = (numel(tf)-1)/2; dk = k;
    tf(dd+1-dk : dd+1+dk) = tf(dd+1-dk : dd+1+dk) + coef * yk;
    if k < p-1, yk = conv(yk, y); end
end
s  = specfact_minphase(tf(:), 2^16);
hb = 1;
for i = 1:p, hb = conv(hb, [1 1]/2); end
h0 = conv(s(:).', hb);
if sum(h0) < 0, h0 = -h0; end
end
