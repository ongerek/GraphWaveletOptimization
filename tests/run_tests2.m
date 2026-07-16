function fails = run_tests2
%RUN_TESTS2 Verification harness for the Wiener-lifting module.
% path bootstrap: locate src/ relative to this file
if exist('design_compaction_fb', 'file') ~= 2
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
end

tr = @(name, err, tol) fprintf('%-52s %.3e  [%s]\n', name, err, ok(err<tol));
fails = 0;

% T1: PR of extracted equivalent filters, random lifting --------------
try, rng(11); catch, rand('seed',11); end %#ok
p = randn(3,1); u = randn(2,1);
kp = -1:1; ku = -1:0;
[h0,h1,g0,g1,off] = lifting_filters(p, kp, u, ku);
e = pr_err(h0,h1,g0,g1,off);
tr('T1 PR of extracted filters (random p,u)', e, 1e-12);
fails = fails + (e>=1e-12);

% T2: LeGall 5/3 recovery: p=[1/2 1/2], u=[1/4 1/4] -------------------
[h0,h1,g0,g1] = lifting_filters([0.5;0.5], 0:1, [0.25;0.25], -1:0);
h0 = h0 / max(abs(h0)) * 6;                     % scale-free compare
e1 = max(abs(h0 - [-1 2 6 2 -1]));
h1 = h1 / max(abs(h1)) * 1;
e2 = max(abs(h1 - [-0.5 1 -0.5]));
tr('T2 LeGall 5/3 h0,h1 recovered', max(e1,e2), 1e-12);
fails = fails + (max(e1,e2)>=1e-12);

% T3: AR(1) Markov blanket: long Wiener predictor collapses -----------
rho = 0.95;
r = acf_from_psd(@(w) psd_ar([1 -rho], 1-rho^2, w));
d8 = design_wiener_lifting(r, 8, 2, struct('dualvm', false));
pth = rho/(1+rho^2);
pref = zeros(8,1); pref(4:5) = pth;             % kp = -3:4 -> k=0,1 at 4,5
e = max(abs(d8.p - pref));
tr('T3 AR(1): 8-tap Wiener predictor = 2-tap Markov', e, 1e-9);
fails = fails + (e>=1e-9);
sdth = (1-rho^2)/(1+rho^2);
e = abs(d8.sigma_d2 - sdth);
tr('T3 AR(1): sigma_d^2 = (1-rho^2)/(1+rho^2)', e, 1e-9);
fails = fails + (e>=1e-9);

% T4: GMRF conditional bound matches AR(1) closed form ----------------
e = abs(gmrf_detail_bound(r) - sdth);
tr('T4 GMRF bound (matrix conditioning) vs closed form', e, 1e-9);
fails = fails + (e>=1e-9);

% T5: Gouze update solves its own quadratic ---------------------------
dw = design_wiener_lifting(r, 4, 4, struct('primalvm', false));
% brute perturbation check: J(u*) <= J(u* + eps e_i)
Jf = @(u) lift_J(u, dw);
J0 = Jf(dw.u); worst = -Inf;
for i = 1:4
    ei = zeros(4,1); ei(i) = 1e-4;
    worst = max(worst, max(J0 - Jf(dw.u+ei), J0 - Jf(dw.u-ei)));
end
tr('T5 Gouze update: local optimality (perturbation)', max(worst,0), 1e-12);
fails = fails + (worst>=1e-12);

% T6: vanishing-moment mechanics --------------------------------------
d22 = design_wiener_lifting(r, 2, 2, struct());
x = ones(64,1);
[~, dd] = lift_apply(x, d22.p, d22.kp, d22.u, d22.ku);
e1 = max(abs(dd));                              % dual VM: d(const)=0
e2 = abs(sum(d22.g1));                          % primal VM: mean(g1)=0
tr('T6 sum(p)=1 -> d(const)=0;  sum(u)=1/2 -> sum(g1)=0', max(e1,e2), 1e-12);
fails = fails + (max(e1,e2)>=1e-12);

% T7: Monte-Carlo validation of analytic sigma_d^2 --------------------
nsim = 2^17;
e = randn(nsim,1); x = filter(1, [1 -rho], e*sqrt(1-rho^2));
[~, dd] = lift_apply(x(1001:end-mod(numel(x)-1000,2)), d22.p, d22.kp, d22.u, d22.ku);
e = abs(var(dd) - d22.sigma_d2)/d22.sigma_d2;
tr('T7 Monte-Carlo sigma_d^2 (rel err, ~1/sqrt(N))', e, 0.03);
fails = fails + (e>=0.03);

fprintf('\n%d test group failures.\n', fails);
end

function s = ok(b), if b, s='PASS'; else, s='FAIL'; end, end

function [a, d] = lift_apply(x, p, kp, u, ku)
x = x(1:2*floor(numel(x)/2));
s = x(1:2:end); t = x(2:2:end);
d = t;
for i = 1:numel(p), d = d - p(i)*circshift(s, -kp(i)); end
a = s;
for j = 1:numel(u), a = a + u(j)*circshift(d, -ku(j)); end
end

function J = lift_J(u, dw)
% analytic Gouze objective rebuilt from stored rd (independent path)
rd = dw.rd; kp = dw.kp; ku = dw.ku; p = dw.p;
np = numel(p); nu = numel(u);
Rdu = rd(abs(ku.' - ku) + 1);
mv  = (min(kp)+min(ku)) : (max(kp)+max(ku));
Cv  = zeros(numel(mv), nu);
for i = 1:np
    for j = 1:nu
        m = kp(i)+ku(j);
        Cv(m-mv(1)+1, j) = Cv(m-mv(1)+1, j) + p(i);
    end
end
Rdv = rd(abs(mv.' - mv) + 1);
rhov = rd(abs(mv.') + 1);
J = u.'*(Rdu + Cv.'*Rdv*Cv)*u - 2*(Cv.'*rhov).'*u + rd(1);
end

function e = pr_err(h0, h1, g0, g1, off)
n = 128; x = randn(n, 1);
a = anab(h0, off(1), x); d = anab(h1, off(2), x);
xh = synb(g0, off(3), a, n) + synb(g1, off(4), d, n);
e = max(abs(xh - x));
end
function a = anab(h, off, x)
n = numel(x); a = zeros(n/2, 1);
for m = 0:n/2-1
    acc = 0;
    for i = 1:numel(h)
        acc = acc + h(i) * x(mod(2*m - (off+i-1), n) + 1);
    end
    a(m+1) = acc;
end
end
function xh = synb(g, off, a, n)
xh = zeros(n, 1);
for m = 0:n/2-1
    for i = 1:numel(g)
        t = mod(2*m + off + i - 1, n) + 1;
        xh(t) = xh(t) + g(i) * a(m+1);
    end
end
end
