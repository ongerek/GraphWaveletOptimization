function out = design_wiener_lifting(r, np, nu, opts)
%DESIGN_WIENER_LIFTING Statistically optimal biorthogonal lifting FB.
%
%   out = DESIGN_WIENER_LIFTING(r, np, nu) designs a predict/update
%   lifting scheme for a WSS input with autocorrelation r
%   (r(1) = r[0]):
%
%     split:    s[n] = x[2n],  t[n] = x[2n+1]
%     predict:  d[n] = t[n] - sum_k p_k s[n+k]        (np taps)
%     update:   a[n] = s[n] + sum_k u_k d[n+k]        (nu taps)
%
%   Perfect reconstruction is STRUCTURAL (any p, u). The design is:
%
%   * p = half-rate WIENER predictor: minimize E[d^2]
%       (R_ee) p = r_eo,  (R_ee)_{kl} = r[2|k-l|],  (r_eo)_k = r[|2k-1|],
%     optionally constrained sum(p) = 1  -> one dual vanishing moment
%     (analysis highpass annihilates constants).
%
%   * u = GOUZE-type update: minimize the reconstruction MSE when the
%     detail branch is discarded (multiresolution approximation quality):
%       J(u) = E[|U d|^2] + E[|d - P U d|^2]
%     which is quadratic in u:  J = u'A u - 2 b'u + r_d[0], with
%       A = R_d + C_v' R_d C_v,   b = C_v' rho_d,   v = p * u (composite),
%     built from the EXACT detail autocorrelation r_d.
%     Optionally constrained sum(u) = 1/2 -> one primal vanishing moment
%     (synthesis wavelet has zero mean; DC routed entirely to a).
%
%   opts fields (optional): .dualvm (default true), .primalvm (default
%   true), .u_fixed (bypass update design, e.g. [1/4 1/4]),
%   .p_fixed (bypass predict design).
%
%   Output struct: p, u, tap offsets kp, ku; equivalent filters
%   h0, h1, g0, g1 with first-tap lags off = [o_h0 o_h1 o_g0 o_g1];
%   sigma_d2 (exact detail variance), sigma_a2, J (update objective).
%
%   See: Gouze, Antonini, Barlaud, Macq, IEEE T-IP 13(12), 2004;
%        Claypoole, Davis, Sweldens, Baraniuk, IEEE T-IP 12(12), 2003.

if nargin < 4, opts = struct(); end
dualvm   = getdef(opts, 'dualvm', true);
primalvm = getdef(opts, 'primalvm', true);
r = r(:);

kp = -floor((np-1)/2) : ceil((np-1)/2);      % predict offsets (even grid)
ku = -ceil(nu/2) : (floor(nu/2) - 1);        % update offsets

% ---------------- predict: constrained Wiener --------------------------
if isfield(opts, 'p_fixed')
    p = opts.p_fixed(:);
else
    Rp = r(2*abs(kp.' - kp) + 1);
    qp = r(abs(2*kp.' - 1) + 1);
    if dualvm
        K  = [Rp, ones(np,1); ones(1,np), 0];
        pz = K \ [qp; 1];
        p  = pz(1:np);
    else
        p  = Rp \ qp;
    end
end

% ---------------- exact detail autocorrelation -------------------------
% d[n] = sum_t w_t x[2n+1-t]:  w_0 = 1,  w_{1-2k} = -p_k
tmin = min(0, 1-2*max(kp)); tmax = max(0, 1-2*min(kp));
w = zeros(1, tmax - tmin + 1);
w(0 - tmin + 1) = 1;
for i = 1:np
    w(1 - 2*kp(i) - tmin + 1) = w(1 - 2*kp(i) - tmin + 1) - p(i);
end
phi = conv(w, fliplr(w));                    % autocorr of w, center = len(w)
cph = numel(w);
Md  = 4*(np + nu) + 8;                       % detail lags needed
rd  = zeros(Md+1, 1);
for m = 0:Md
    j  = -(numel(w)-1) : (numel(w)-1);
    rd(m+1) = phi(cph + j) * r(abs(2*m - j) + 1);
end
sigma_d2 = rd(1);

% ---------------- update: Gouze quadratic -------------------------------
if isfield(opts, 'u_fixed')
    u = opts.u_fixed(:); J = NaN;
else
    Rdu = rd(abs(ku.' - ku) + 1);
    % composite lags v = p * u:  v_m, m = kp_i + ku_j
    mv  = (min(kp)+min(ku)) : (max(kp)+max(ku));
    Cv  = zeros(numel(mv), nu);
    for i = 1:np
        for j = 1:nu
            m = kp(i) + ku(j);
            Cv(m - mv(1) + 1, j) = Cv(m - mv(1) + 1, j) + p(i);
        end
    end
    Rdv  = rd(abs(mv.' - mv) + 1);
    rhov = rd(abs(mv.') + 1);
    A = Rdu + Cv.'*Rdv*Cv;
    b = Cv.'*rhov;
    if primalvm
        K  = [A, ones(nu,1); ones(1,nu), 0];
        uz = K \ [b; 0.5];
        u  = uz(1:nu);
    else
        u  = A \ b;
    end
    J = u.'*A*u - 2*b.'*u + rd(1);
end

% ---------------- equivalent filters (numeric extraction) --------------
[h0, h1, g0, g1, off] = lifting_filters(p, kp, u, ku);

% exact subband variances from equivalent analysis filters
sigma_a2 = eqvar(h0, r);

out = struct('p', p, 'u', u, 'kp', kp, 'ku', ku, ...
             'h0', h0, 'h1', h1, 'g0', g0, 'g1', g1, 'off', off, ...
             'sigma_d2', sigma_d2, 'sigma_a2', sigma_a2, 'J', J, ...
             'rd', rd);
end

function v = getdef(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end

function s2 = eqvar(h, r)
n = numel(h); q = conv(h(:).', fliplr(h(:).')); c = n;
s2 = q(c)*r(1) + 2*(q(c+1:end) * r(2:n));
end
