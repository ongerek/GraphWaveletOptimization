function out = design_compaction_fb(r, N, p, opts)
%DESIGN_COMPACTION_FB Signal-adapted orthonormal two-channel FIR filterbank.
%
%   out = DESIGN_COMPACTION_FB(r, N, p) designs a length-N (N even)
%   real orthonormal (paraunitary/CQF) lowpass filter h0 with p vanishing
%   moments that MAXIMIZES the lowpass subband variance
%
%       sigma_a^2 = sum_k f[k] r[k],      f = h0 * h0(-.)   (product filter)
%
%   for a zero-mean WSS input with autocorrelation sequence r
%   (r(1) = r[0], r(2) = r[1], ...; need length(r) >= 2N).
%
%   Method (Kirac & Vaidyanathan 1998 / Moulin et al. 1997, LP variant):
%     F(z) = H0(z) H0(1/z) is parametrized as F = B(z)^p T(z), where
%     B(z) = (1+z)(1+z^{-1})/4  (so |H0|^2 has a zero of order 2p at pi),
%     T(z) symmetric of degree d = N-1-p.  Then:
%       * orthonormality (half-band: f[0]=1, f[2k]=0)  -> LINEAR equalities
%       * F >= 0  <=>  T(w) >= 0                        -> LINEAR inequalities
%                                                          (discretized grid)
%       * sigma_a^2                                     -> LINEAR objective
%     i.e. a linear program (the exact continuum version is an SDP via the
%     trace parametrization; the dense-grid LP + a-posteriori verification
%     is solver-free and accurate to ~1e-9 here).
%     The optimal T is blended with a strictly positive feasible reference
%     (Daubechies-p) with weight opts.blend to guarantee a numerically
%     clean minimum-phase spectral factorization (cepstral method).
%
%   opts fields (all optional):
%     .ngrid  grid size for T(w)>=0 on [0,pi]   (default 4096)
%     .blend  blend weight lambda toward Daubechies-p (default 1e-8)
%     .nfft   FFT size for spectral factorization (default 2^16)
%
%   Output struct:
%     .h0, .h1      analysis filters (orthonormal pair, h1 = CQF of h0)
%     .f            product filter coefficients (length 2N-1, center N)
%     .sigma_a2     achieved lowpass subband variance
%     .sigma_a2_lp  LP optimal value (before blending; equal to ~1e-8)
%     .err_hb       max |half-band violation| of conv(h0,fliplr(h0))
%     .err_orth     max |<h0, h0(.-2k)> - delta[k]|
%     .Tmin         min of T(w) on a fine verification grid
%
%   Requires: LP_SOLVE (wrapper around linprog or GLPK), SPECFACT_MINPHASE.

if nargin < 4, opts = struct(); end
Ng   = getfielddef(opts, 'ngrid', 4096);
lam  = getfielddef(opts, 'blend', 1e-8);
Kfft = getfielddef(opts, 'nfft',  2^16);

r = r(:);
if mod(N,2) ~= 0, error('N must be even.'); end
M = N/2;
if p < 0 || p > M, error('Need 0 <= p <= N/2.'); end
if numel(r) < 2*N, error('Provide autocorrelation to lag >= 2N-1.'); end
d = 2*M - 1 - p;                       % degree of T

% ---- B(z)^p coefficients (centered, length 2p+1) --------------------
b = 1;
for i = 1:p, b = conv(b, [1 2 1]/4); end

% ---- linear maps: t (d+1 free params) -> tf (Laurent) -> f ----------
E = zeros(2*d+1, d+1);                 % symmetric embedding
E(d+1, 1) = 1;
for m = 1:d
    E(d+1-m, m+1) = 1;
    E(d+1+m, m+1) = 1;
end
Cb = convmtx_local(b(:), 2*d+1);       % (4M-1) x (2d+1) convolution matrix
Ct = Cb * E;                           % f = Ct * t,  length Lf = 4M-1
Lf = 4*M - 1;
c0 = 2*M;                              % center index of f (1-based)

% ---- LP data ---------------------------------------------------------
rw   = r(abs((1:Lf).' - c0) + 1);      % r[|k - c0|]
cvec = Ct.' * rw;                      % sigma_a^2 = cvec' * t
Aeq  = Ct(c0 : 2 : c0 + 2*(M-1), :);   % f[0]=1, f[2k]=0 (k=1..M-1)
beq  = [1; zeros(M-1, 1)];
w    = linspace(0, pi, Ng).';
Aineq = -[ones(Ng,1), 2*cos(w * (1:d))];   % -T(w) <= 0
bineq = zeros(Ng, 1);

[t, ok] = lp_solve(-cvec, Aineq, bineq, Aeq, beq);
if ~ok, error('LP solver failed.'); end
sigma_a2_lp = cvec.' * t;

% ---- blend with strictly positive reference (keeps all constraints) --
% T_ref >= 2 on [0,pi]; choose lambda adaptively so that the blended T is
% strictly positive despite LP-grid discretization / solver tolerances.
tref = reference_T(p, d);
wv   = linspace(0, pi, 16*Ng).';
Cosv = [ones(numel(wv),1), 2*cos(wv*(1:d))];
Tmin0 = min(Cosv * t);
lam  = max(lam, 1.5 * max(0, -Tmin0) / 2);
t    = (1 - lam)*t + lam*tref;
tf   = E * t;

% ---- spectral factorization and filter assembly ----------------------
s  = specfact_minphase(tf, Kfft);      % min-phase, |S|^2 = T
hb = 1;
for i = 1:p, hb = conv(hb, [1 1]/2); end
h0 = conv(s(:).', hb);
if sum(h0) < 0, h0 = -h0; end
h1 = fliplr(h0) .* (-1).^(0:N-1);      % CQF highpass

% ---- verification -----------------------------------------------------
f  = conv(h0, fliplr(h0));
ev = f(c0+2 : 2 : end);
err_hb   = max(abs([f(c0) - 1, ev]));
sigma_a2 = f * r(abs((1:Lf).' - c0) + 1);
Tmin = min(Cosv * t);

out = struct('h0', h0, 'h1', h1, 'f', f, ...
             'sigma_a2', sigma_a2, 'sigma_a2_lp', sigma_a2_lp, ...
             'err_hb', err_hb, 'err_orth', err_hb, 'Tmin', Tmin, ...
             't', t, 'p', p, 'N', N, 'lambda', lam);
end

% ----------------------------------------------------------------------
function v = getfielddef(s, name, def)
if isfield(s, name), v = s.(name); else, v = def; end
end

function C = convmtx_local(b, n)
% Convolution matrix: C*x = conv(b, x), b column of length nb.
nb = numel(b);
C  = zeros(nb + n - 1, n);
for j = 1:n
    C(j : j+nb-1, j) = b;
end
end

function t = reference_T(p, d)
% Strictly positive feasible T: Daubechies-p residual polynomial
% T_db(z) = 2 * sum_{k=0}^{p-1} C(p-1+k, k) * (1 - B(z))^k,  degree p-1,
% zero-padded to degree d.  For p = 0: lazy filterbank, T = 1.
if p == 0
    t = [1; zeros(d, 1)];
    return;
end
y  = [-1 2 -1]/4;                      % 1 - B(z), centered
tf = zeros(1, 2*(p-1)+1);
yk = 1;                                % y^0
for k = 0:p-1
    coef = 2 * nchoosek(p-1+k, k);
    tf   = addcentered(tf, coef * yk);
    if k < p-1, yk = conv(yk, y); end
end
dd = p - 1;                            % degree of tf
t  = zeros(d+1, 1);
t(1 : dd+1) = tf(dd+1 : end).';        % nonnegative-lag part
end

function a = addcentered(a, b)
% Add two centered symmetric Laurent coefficient vectors (odd lengths).
da = (numel(a)-1)/2; db = (numel(b)-1)/2;
if db > da, tmp = a; a = b; b = tmp; da2 = da; da = db; db = da2; end
c0 = da + 1;
a(c0-db : c0+db) = a(c0-db : c0+db) + b;
end
