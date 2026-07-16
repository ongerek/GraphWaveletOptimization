function [GdB, info] = coding_gain_bior_adapt(designer, r, L)
%CODING_GAIN_BIOR_ADAPT L-level dyadic coding gain for a (possibly
%   per-level re-designed) biorthogonal FB, Katto-Yasuda weighting.
%   designer: @(r) -> struct with fields h0, h1, g0, g1; called at every
%   level with the autocorrelation of the current lowpass signal.
r = r(:); sigx2 = r(1);
sd2 = zeros(L,1); wd = zeros(L,1);
chain = 1;
for j = 1:L
    fb = designer(r);
    h0 = fb.h0(:).'; h1 = fb.h1(:).';
    q0 = conv(h0, fliplr(h0));  n0 = numel(h0);
    q1 = conv(h1, fliplr(h1));  n1 = numel(h1);
    sd2(j) = subvar(q1, n1, r);
    geq    = conv(chain, upsamp(fb.g1(:).', 2^(j-1)));
    wd(j)  = sum(geq.^2);
    r      = nextacf(q0, n0, r);
    chain  = conv(chain, upsamp(fb.g0(:).', 2^(j-1)));
end
saL2 = r(1);
wa   = sum(chain.^2);
den  = (saL2*wa)^(2^-L) * prod((sd2 .* wd) .^ (2.^-(1:L)'));
GdB  = 10*log10(sigx2 / den);
info = struct('sd2', sd2, 'wd', wd, 'saL2', saL2, 'wa', wa);
end

function v = subvar(q, n, r)
c0 = n;
v  = q(c0)*r(1) + 2*(q(c0+1:c0+n-1) * r(2:n));
end

function rn = nextacf(q, n, r)
Lr = floor((numel(r) - n)/2);
ks = -(n-1):(n-1);
qv = q(n + ks);
rn = zeros(Lr+1,1);
for m = 0:Lr
    rn(m+1) = qv * r(abs(2*m - ks) + 1);
end
end

function y = upsamp(x, f)
if f == 1, y = x; return; end
y = zeros(1, f*(numel(x)-1)+1);
y(1:f:end) = x;
end
