function fails = run_tests
%RUN_TESTS Verification harness for the compaction-fb package.
tol_report = @(name, err, tol) fprintf('%-46s %.3e  [%s]\n', name, err, ...
    merge_str(err < tol));
fails = 0;

% ---- T1: spectral factorization on random positive T ----------------
rng_seed(7);
d = 11;
s_true = randn(1, d+1);
tf = conv(s_true, fliplr(s_true));   % T = S S~, guaranteed >= 0
% make strictly positive by adding margin
tf(d+1) = tf(d+1) + 0.1*sum(abs(tf));
s = specfact_minphase(tf(:), 2^16);
err = max(abs(conv(s(:).', fliplr(s(:).')) - tf));
tol_report('T1 specfact: |S*S~ - T|_inf', err, 1e-10);
fails = fails + (err >= 1e-10);

% ---- T2: Daubechies generator vs known db2 coefficients -------------
h = db_filter(2);
href = [0.482962913144534 0.836516303737808 0.224143868042013 -0.129409522551260];
err = min(max(abs(h - href)), max(abs(h + fliplr(href))));
tol_report('T2 db2 vs literature coefficients', err, 1e-9);
fails = fails + (err >= 1e-9);
% orthonormality of db8
h = db_filter(8);
q = conv(h, fliplr(h)); c0 = numel(h);
err = max(abs([q(c0)-1, q(c0+2:2:end)]));
tol_report('T2 db8 orthonormality (half-band err)', err, 1e-9);
fails = fails + (err >= 1e-9);

% ---- T3: CDF 9/7 -- lengths, PR, literature values -------------------
[h0, h1, g0, g1] = cdf97_filters();
fprintf('%-46s %d/%d/%d/%d\n', 'T3 9/7 lengths (expect 9/7/7/9)', ...
        numel(h0), numel(h1), numel(g0), numel(g1));
h0ref = sqrt(2)*[0.026748757410810 -0.016864118442875 -0.078223266528990 ...
                 0.266864118442875  0.602949018236360  0.266864118442875 ...
                -0.078223266528990 -0.016864118442875  0.026748757410810];
err = max(abs(h0 - h0ref));
tol_report('T3 9/7 h0 vs literature (sqrt2 norm)', err, 1e-9);
fails = fails + (err >= 1e-9);
% PR on random signal via lifting round trip is structural; check the
% extracted FILTERS satisfy PR:  sum_n g0[n-2k] h0[n] + g1[n-2k] h1[n] ...
% direct operator test on a random vector with the extracted filters:
[h0, h1, g0, g1, off] = cdf97_filters();
err = pr_error(h0, h1, g0, g1, off);
tol_report('T3 9/7 PR (polyphase determinant test)', err, 1e-10);
fails = fails + (err >= 1e-10);

% ---- T4: LP design on AR(1): must match/beat Daubechies -------------
rho = 0.95;
r = acf_from_psd(@(w) psd_ar([1 -rho], 1-rho^2, w));
res = design_compaction_fb(r, 16, 2);
tol_report('T4 adapted N=16 half-band error', res.err_hb, 1e-8);
fails = fails + (res.err_hb >= 1e-8);
fprintf('%-46s %.6f (Tmin = %.2e)\n', 'T4 sigma_a^2 (LP/achieved)', ...
        res.sigma_a2, res.Tmin);
h8 = db_filter(8);
f8 = conv(h8, fliplr(h8)); c8 = numel(h8);
sa_db = f8(c8)*r(1) + 2*(f8(c8+1:c8+15) * r(2:16));
fprintf('%-46s %.6f vs db8 %.6f  [%s]\n', 'T4 compaction: adapted >= db8', ...
        res.sigma_a2, sa_db, merge_str(res.sigma_a2 >= sa_db - 1e-9));
fails = fails + (res.sigma_a2 < sa_db - 1e-9);

% ---- T5: coding gain sanity, AR(1) rho=0.95 --------------------------
G = coding_gain_ortho(@(rr) h8, r, 5);
bound = -10*log10(1 - rho^2);   % distortion-rate bound: 1/(1-rho^2)
fprintf('%-46s %.3f dB (bound %.3f dB) [%s]\n', ...
        'T5 db8 5-level gain vs AR(1) bound', G, bound, ...
        merge_str(G < bound && G > 9));
fails = fails + ~(G < bound && G > 9);

% ---- T6: bior gain formula reduces to ortho for orthonormal FB ------
h1o = fliplr(h8).*(-1).^(0:15);
Gb = coding_gain_bior(h8, h1o, fliplr(h8), fliplr(h1o), r, 4);
Go = coding_gain_ortho(@(rr) h8, r, 4);
tol_report('T6 bior formula == ortho formula (dB)', abs(Gb-Go), 1e-9);
fails = fails + (abs(Gb-Go) >= 1e-9);

fprintf('\n%d test group failures.\n', fails);

% ======================================================================
function s = merge_str(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
function rng_seed(k)
try, rng(k); catch, rand('seed', k); randn('seed', k); end %#ok
end
function e = pr_error(h0, h1, g0, g1, off)
% Exact time-domain PR test on a random periodic signal using the
% extracted filters WITH their tap offsets:
%   a[m] = sum_k h0a[k] x[2m-k],  xhat[t] = sum_m g0a[t-2m] a[m] + ...
n = 128; x = randn(n, 1);
a = anab(h0, off(1), x);  d = anab(h1, off(2), x);
xh = synb(g0, off(3), a, n) + synb(g1, off(4), d, n);
e = max(abs(xh - x));
end
function a = anab(h, off, x)
n = numel(x); a = zeros(n/2, 1);
for m = 0:n/2-1
    acc = 0;
    for i = 1:numel(h)
        k = off + i - 1;
        acc = acc + h(i) * x(mod(2*m - k, n) + 1);
    end
    a(m+1) = acc;
end
end
function xh = synb(g, off, a, n)
xh = zeros(n, 1);
for m = 0:n/2-1
    for i = 1:numel(g)
        k = off + i - 1;
        t = mod(2*m + k, n) + 1;
        xh(t) = xh(t) + g(i) * a(m+1);
    end
end
end

end
