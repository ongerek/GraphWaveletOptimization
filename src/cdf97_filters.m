function [h0, h1, g0, g1, off] = cdf97_filters()
%CDF97_FILTERS CDF 9/7 biorthogonal analysis/synthesis filters, built
%   from the lifting factorization (Daubechies-Sweldens), so perfect
%   reconstruction is guaranteed by construction. Normalized so that
%   sum(h0) = sqrt(2) and H1(pi) = -sqrt(2) (textbook convention).
%   off = [off_h0 off_h1 off_g0 off_g1]: lag of the first tap of each
%   filter, so that  a[m] = sum_k h0(k-off+1) x[2m - k]  etc., and PR is
%   x[t] = sum_m g0(.) a[m] + g1(.) d[m] with the same lag convention.
n  = 64;                               % periodic extraction length
m0 = n/4;                              % probe output index (0-based)
% ---- analysis filters: impulse at time t -> a[m] = h0[2m - t] --------
h0 = zeros(1, 2*n); h1 = zeros(1, 2*n);   % index k = -n..n-1 -> k+n+1
for t = 0:n-1
    x = zeros(n, 1); x(t+1) = 1;
    [a, d] = lift97_fwd(x);
    k = 2*m0 - t;
    if abs(k) <= n-1
        h0(k+n+1) = a(m0+1);
        h1(k+n+1) = d(m0+1);
    end
end
% ---- synthesis: a = delta_{m0}, d = 0  ->  x[t] = g0[t - 2 m0] -------
av = zeros(n/2, 1); dv = zeros(n/2, 1);
av(m0+1) = 1;
x = lift97_inv(av, dv);
g0 = readout(x, n, m0);
av(:) = 0; dv(m0+1) = 1;
x = lift97_inv(av, dv);
g1 = readout(x, n, m0);
[h0, o0] = trim(h0, n); [h1, o1] = trim(h1, n);
[g0, o2] = trim(g0, n); [g1, o3] = trim(g1, n);
off = [o0 o1 o2 o3];        % lag index of first tap of each filter
% ---- normalization (coding gain is invariant to this rescaling) ------
s0 = sqrt(2) / sum(h0);
h0 = s0 * h0;  g0 = g0 / s0;
s1 = -sqrt(2) / (((-1).^(0:numel(h1)-1)) * h1(:));
h1 = s1 * h1;  g1 = g1 / s1;
end

function [s, t] = lift97_fwd(x)
[al, be, ga, de] = coeffs97();
s = x(1:2:end); t = x(2:2:end);
t = t + al*(s + circshift(s, -1));
s = s + be*(t + circshift(t,  1));
t = t + ga*(s + circshift(s, -1));
s = s + de*(t + circshift(t,  1));
end

function x = lift97_inv(s, t)
[al, be, ga, de] = coeffs97();
s = s - de*(t + circshift(t,  1));
t = t - ga*(s + circshift(s, -1));
s = s - be*(t + circshift(t,  1));
t = t - al*(s + circshift(s, -1));
x = zeros(2*numel(s), 1);
x(1:2:end) = s; x(2:2:end) = t;
end

function [al, be, ga, de] = coeffs97()
al = -1.586134342059924;  be = -0.052980118572961;
ga =  0.882911075530934;  de =  0.443506852043971;
end

function g = readout(x, n, m0)
g = zeros(1, 2*n);
for t = 0:n-1
    k = t - 2*m0;
    if abs(k) <= n-1, g(k+n+1) = x(t+1); end
end
end

function [h, off] = trim(h, n)
i   = find(abs(h) > 1e-12);
off = i(1) - (n + 1);       % first-tap lag (index k stored at k+n+1)
h   = h(i(1):i(end));
end
