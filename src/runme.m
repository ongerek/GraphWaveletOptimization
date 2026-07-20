function runme
%RUNME Regenerate every insertable asset of the manuscript from the
%   GraphWaveletOptimization package: Figures 1-4 (PNG 300 dpi + EPS) and the data rows of
%   Tables 1-7 as ready-to-paste LaTeX fragments.
%
%   Usage: place this file in the GraphWaveletOptimization package folder (with
%   design_compaction_fb.m etc. on the path) and run
%       >> runme
%   Outputs are written to ./paper_out/ :
%     fig1_bimodal_responses.{png,eps}   -> manuscript Fig. 1
%     fig2_predictor_taps.{png,eps}      -> manuscript Fig. 2
%     fig3_graph_pyramid.{png,eps}       -> manuscript Fig. 3
%     fig4_localization_sweep.{png,eps}  -> manuscript Fig. 8 (k-sweep)
%     table1_gains_module1.tex           -> Table 1 body rows
%     table2_compaction.tex              -> Table 2 body rows
%     table3_gains_module2.tex           -> Table 3 body rows
%     table4_vm_price.tex                -> Table 4 body rows
%     table5_path.tex                    -> Table 5 body rows
%     table6_sensor.tex                  -> Table 6 body rows
%     table7_learning.tex                -> Table 7 body rows
%   All numbers are also echoed to the console. Figure blocks are
%   wrapped in try/catch so the numeric outputs complete even on a
%   display-less machine. Runtime: a few minutes.

od = 'paper_out';
if ~exist(od, 'dir'), mkdir(od); end
seed(0);

% ---------------------------------------------------------------------
% Processes P1-P3 (unit variance), L = 4 levels
% ---------------------------------------------------------------------
L = 4; K = 2^15; N = 16;
S1 = @(w) psd_ar([1 -0.95], 1 - 0.95^2, w);
rho2 = 0.975; th2 = 0.8*pi;
S2 = @(w) psd_ar([1, -2*rho2*cos(th2), rho2^2], 1, w);
bump = @(w,w0,s) exp(-0.5*((w-w0)/s).^2) + exp(-0.5*((2*pi-w-w0)/s).^2);
S3 = @(w) bump(w,0.15*pi,0.06*pi) + 0.8*bump(w,0.70*pi,0.05*pi) + 0.01;
Sfun  = {S1, S2, S3};
pname = {'P1 AR(1) $\rho{=}0.95$', 'P2 AR(2) $@0.8\pi$', 'P3 bimodal'};
w = 2*pi*(0:K-1).'/K;

hdb = {db_filter(2), db_filter(4), db_filter(8)};
[b0, b1, bg0, bg1] = cdf97_filters();
h8  = hdb{3};
h8s = struct('h0', h8, 'h1', fliplr(h8).*(-1).^(0:15), ...
             'g0', fliplr(h8), 'g1', (-1).^(0:15).*h8);
[l0h, l1h, lg0, lg1] = lifting_filters([.5;.5], 0:1, [.25;.25], -1:0);
l53 = struct('h0', l0h, 'h1', l1h, 'g0', lg0, 'g1', lg1);

% =====================================================================
% TABLES 1 & 2  (Module 1)
% =====================================================================
fprintf('\n== Tables 1 & 2 (Module 1) ==\n');
f1 = fopen(fullfile(od, 'table1_gains_module1.tex'), 'w');
f2 = fopen(fullfile(od, 'table2_compaction.tex'), 'w');
for ip = 1:3
    r = acf_from_psd(Sfun{ip}, K); r = r / r(1);
    G = zeros(1, 7);
    for i = 1:3, G(i) = coding_gain_ortho(@(rr) hdb{i}, r, L); end
    G(4) = coding_gain_bior(b0, b1, bg0, bg1, r, L);
    des  = design_compaction_fb(r, N, 2);
    G(5) = coding_gain_ortho(@(rr) des.h0, r, L);
    G(6) = coding_gain_ortho(@(rr) adap_orth(rr, N, 2), r, L);
    G(7) = coding_gain_ortho(@(rr) adap_orth(rr, N, 0), r, L);
    Sg = Sfun{ip}(w); Sg = Sg / mean(Sg);
    Gmax = -10*log10(exp(mean(log(Sg))));
    fprintf(f1, '%s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n', ...
            pname{ip}, G, Gmax);
    % single-stage compaction (Table 2)
    fq = conv(h8, fliplr(h8)); c8 = numel(h8);
    sa_db = fq(c8)*r(1) + 2*(fq(c8+1:end) * r(2:c8));
    wf = linspace(0, pi, 2^14).';
    Sw = Sfun{ip}(wf); Sm = Sfun{ip}(pi - wf);
    sa_id = trapz(wf, max(Sw, Sm)) / trapz(wf, Sw);
    fprintf(f2, '%s & %.4f & %.4f & %.4f\\\\\n', ...
            pname{ip}, sa_db, des.sigma_a2, sa_id);
    fprintf('  %-24s gains done (bound %.3f dB)\n', pname{ip}, Gmax);
end
fclose(f1); fclose(f2);

% =====================================================================
% FIGURE 1  (bimodal responses)
% =====================================================================
try
    r  = acf_from_psd(S3, K); r = r/r(1);
    d2 = design_compaction_fb(r, N, 2);
    d0 = design_compaction_fb(r, N, 0);
    wf = linspace(0, pi, 1024).';
    Hm = @(h) abs(exp(-1i*wf*(0:numel(h)-1)) * h(:)).^2;
    Sw = S3(wf); Sw = Sw / max(Sw) * 2;
    fh = figure;
    plot(wf/pi, Sw, 'k-', 'linewidth', 1.6); hold on;
    plot(wf/pi, Hm(h8),    'b--', 'linewidth', 1.3);
    plot(wf/pi, Hm(d2.h0), 'r-',  'linewidth', 1.5);
    plot(wf/pi, Hm(d0.h0), 'm-.', 'linewidth', 1.3);
    grid on; xlabel('\omega / \pi'); ylim([-0.05 2.25]);
    legend('PSD (scaled)', 'db8 |H_0|^2', 'adapted p=2 |H_0|^2', ...
           'adapted p=0 |H_0|^2', 'location', 'north');
    saveboth(fh, od, 'fig1_bimodal_responses');
catch err
    warning('Fig 1 skipped: %s', err.message);
end

% =====================================================================
% TABLE 3  (Module 2 gains)
% =====================================================================
fprintf('\n== Table 3 (Module 2) ==\n');
f3 = fopen(fullfile(od, 'table3_gains_module2.tex'), 'w');
ofree = struct('dualvm', false, 'primalvm', false);
for ip = 1:3
    r = acf_from_psd(Sfun{ip}, K); r = r/r(1);
    G = zeros(1, 8);
    G(1) = coding_gain_bior_adapt(@(rr) l53, r, L);
    G(2) = coding_gain_bior(b0, b1, bg0, bg1, r, L);
    G(3) = coding_gain_bior_adapt(@(rr) h8s, r, L);
    G(4) = coding_gain_ortho(@(rr) adap_orth(rr, N, 2), r, L);
    w22  = design_wiener_lifting(r, 2, 2);
    G(5) = coding_gain_bior_adapt(@(rr) fbw(w22), r, L);
    w66  = design_wiener_lifting(r, 6, 6);
    G(6) = coding_gain_bior_adapt(@(rr) fbw(w66), r, L);
    G(7) = coding_gain_bior_adapt(@(rr) fbw(design_wiener_lifting(rr,6,6)), r, L);
    G(8) = coding_gain_bior_adapt(@(rr) fbw(design_wiener_lifting(rr,20,6,ofree)), r, L);
    fprintf(f3, '%s & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n', ...
            pname{ip}, G);
    fprintf('  %-24s done\n', pname{ip});
end
fclose(f3);

% =====================================================================
% TABLE 4 + FIGURE 2  (VM price on P2; predictor taps)
% =====================================================================
fprintf('\n== Table 4 / Fig 2 (VM price on P2) ==\n');
r2 = acf_from_psd(S2, K); r2 = r2 / r2(1);
nps = [2 6 12 20 32];
sc = zeros(size(nps)); sn = sc;
for i = 1:numel(nps)
    dc = design_wiener_lifting(r2, nps(i), 4);
    dn = design_wiener_lifting(r2, nps(i), 4, struct('dualvm', false));
    sc(i) = dc.sigma_d2; sn(i) = dn.sigma_d2;
end
bnd = gmrf_detail_bound(r2);
f4 = fopen(fullfile(od, 'table4_vm_price.tex'), 'w');
fprintf(f4, 'constrained   & %.5f & %.5f & %.5f & %.5f & %.5f\\\\\n', sc);
fprintf(f4, 'unconstrained & %.5f & %.5f & %.5f & %.5f & %.5f\\\\\n', sn);
fprintf(f4, '%% GMRF bound (quote in caption): %.5f\n', bnd);
fclose(f4);
fprintf('  bound %.5f | constrained(32) %.5f | unconstrained(6) %.5f\n', ...
        bnd, sc(end), sn(2));
try
    dV = design_wiener_lifting(r2, 20, 4);
    dF = design_wiener_lifting(r2, 20, 4, struct('dualvm', false));
    fh = figure;
    stem(dV.kp, dV.p, 'filled'); hold on;
    stem(dF.kp, dF.p, 'r');
    grid on; xlabel('tap offset k (even-grid)'); ylabel('p_k');
    legend('constrained  \Sigma p_k = 1', 'unconstrained', ...
           'location', 'northeast');
    saveboth(fh, od, 'fig2_predictor_taps');
catch err
    warning('Fig 2 skipped: %s', err.message);
end

% =====================================================================
% GRAPH SETUP (path + sensor, identical to demo_graph_benchmark)
% =====================================================================
fprintf('\n== Tables 5 & 6 / Figs 3 & 4 (Module 3) ==\n');
rho = 0.95; n = 256;
Jp = (diag([1, (1+rho^2)*ones(1,n-2), 1]) + ...
      diag(-rho*ones(1,n-1),1) + diag(-rho*ones(1,n-1),-1)) / (1-rho^2);
Sp = inv(Jp);
Lp = diag([1, 2*ones(1,n-2), 1]) - diag(ones(1,n-1),1) - diag(ones(1,n-1),-1);
Gs = sensor_graph_demo(200, 11);
Ss = inv(Gs.L + 0.1*eye(200)); Ss = Ss / mean(diag(Ss));

% ---- Table 5: path ----------------------------------------------------
f5 = fopen(fullfile(od, 'table5_path.tex'), 'w');
pts = {'average', 'harmonic', 'wiener', 'exact'};
lbl = {'graph lifting, predict = average', ...
       'graph lifting, predict = harmonic', ...
       'graph lifting, predict = Wiener 1-hop', ...
       'graph lifting, predict = exact'};
for i = 1:4
    pyr = graph_lifting_pyramid(Lp, Sp, 4, ...
          struct('predict', pts{i}, 'hops', 1, 'update', 'gouze'));
    g = tc_gain(pyr.T, pyr.S, Sp);
    fprintf(f5, '%s & %.3f\\\\\n', lbl{i}, g);
    fprintf('  path %-10s %.3f dB\n', pts{i}, g);
end
rp = acf_from_psd(S1, K); rp = rp / rp(1);
Gc = coding_gain_bior_adapt(@(rr) fbw(design_wiener_lifting(rp, 2, 2)), rp, 4);
fprintf(f5, 'classical Module-2 WL(2,2) (stationary, no boundary) & %.3f\\\\\n', Gc);
lam = eig(Sp);
fprintf(f5, 'GFT = KLT bound (AM/GM of eigenvalues) & %.3f\\\\\n', ...
        10*log10(mean(lam)/exp(mean(log(lam)))));
fclose(f5);

% ---- Table 6: sensor network -----------------------------------------
f6 = fopen(fullfile(od, 'table6_sensor.tex'), 'w');
cfg = {'average', 0, 'average (topology only)';
       'harmonic', 0, 'harmonic';
       'wiener',   1, 'Wiener 1-hop';
       'wiener',   2, 'Wiener 2-hop';
       'exact',    0, 'exact conditional mean'};
for i = 1:size(cfg, 1)
    g = zeros(1, 2); ut = {'none', 'gouze'};
    for k = 1:2
        pyr = graph_lifting_pyramid(Gs.L, Ss, 3, struct('predict', cfg{i,1}, ...
              'hops', max(cfg{i,2},1), 'update', ut{k}));
        g(k) = tc_gain(pyr.T, pyr.S, Ss);
    end
    fprintf(f6, '%s & %.3f & %.3f\\\\\n', cfg{i,3}, g);
    fprintf('  sensor %-24s %.3f | %.3f dB\n', cfg{i,3}, g);
end
pyl = graph_lifting_pyramid(Gs.L, Ss, 3, struct('predict','wiener', ...
      'hops',1,'update','gouze','uhops',1));
gl = tc_gain(pyl.T, pyl.S, Ss);
fprintf(f6, 'Wiener 1-hop + 1-hop-masked Gouze & \\multicolumn{2}{c}{%.3f}\\\\\n', gl);
lam = eig(Ss);
fprintf(f6, 'GFT = KLT bound & \\multicolumn{2}{c}{%.3f}\\\\\n', ...
        10*log10(mean(lam)/exp(mean(log(lam)))));
fclose(f6);

% ---- Figure 3: graph, bipartition, Kron levels ------------------------
try
    e1 = graph_bipartition(Gs.L);
    L1 = kron_reduce(Gs.L, e1, 1e-3);  xy1 = Gs.xy(e1, :);
    e2 = graph_bipartition(L1);
    L2 = kron_reduce(L1, e2, 1e-3);    xy2 = xy1(e2, :);
    fh = figure('position', [100 100 1200 380]);
    subplot(1,3,1);
    gplot(Gs.A, Gs.xy, '-'); hold on;
    plot(Gs.xy(e1,1),  Gs.xy(e1,2),  'bo', 'markerfacecolor', 'b', 'markersize', 4);
    plot(Gs.xy(~e1,1), Gs.xy(~e1,2), 'rs', 'markerfacecolor', 'r', 'markersize', 4);
    axis square off; title('level 0 + bipartition (kept = blue)');
    subplot(1,3,2);
    gplot(adjof(L1), xy1, '-'); hold on;
    plot(xy1(:,1), xy1(:,2), 'bo', 'markerfacecolor', 'b', 'markersize', 4);
    axis square off; title(sprintf('Kron level 1 (n = %d)', size(L1,1)));
    subplot(1,3,3);
    gplot(adjof(L2), xy2, '-'); hold on;
    plot(xy2(:,1), xy2(:,2), 'bo', 'markerfacecolor', 'b', 'markersize', 4);
    axis square off; title(sprintf('Kron level 2 (n = %d)', size(L2,1)));
    saveboth(fh, od, 'fig3_graph_pyramid');
catch err
    warning('Fig 3 skipped: %s', err.message);
end

% ---- Figure 4: localization sweep (manuscript Fig. 8) -----------------
try
    kk = [1 2 3 4 6 8 10 14 18 24];
    g0 = nan(size(kk)); gg = nan(size(kk));
    for i = 1:numel(kk)
        p0 = graph_lifting_pyramid(Gs.L, Ss, 3, struct('predict','wiener_knn', ...
             'knn', kk(i), 'update', 'none'));
        pg = graph_lifting_pyramid(Gs.L, Ss, 3, struct('predict','wiener_knn', ...
             'knn', kk(i), 'update', 'gouze'));
        g0(i) = tc_gain(p0.T, p0.S, Ss);
        gg(i) = tc_gain(pg.T, pg.S, Ss);
    end
    px0 = graph_lifting_pyramid(Gs.L, Ss, 3, struct('predict','exact','update','none'));
    pxg = graph_lifting_pyramid(Gs.L, Ss, 3, struct('predict','exact','update','gouze'));
    gx0 = tc_gain(px0.T, px0.S, Ss);  gxg = tc_gain(pxg.T, pxg.S, Ss);
    e1 = graph_bipartition(Gs.L);
    Ad = double(Gs.W > 0); Ad(1:size(Ad,1)+1:end) = 0;
    R = Ad; kh = zeros(1,2); gh0 = zeros(1,2);
    for h = 1:2
        if h == 2, R = double((R*Ad + R) > 0); end
        kh(h) = mean(sum(R(~e1, e1) > 0, 2));
        ph = graph_lifting_pyramid(Gs.L, Ss, 3, struct('predict','wiener', ...
             'hops', h, 'update', 'none'));
        gh0(h) = tc_gain(ph.T, ph.S, Ss);
    end
    lam = eig(Ss); gk = 10*log10(mean(lam)/exp(mean(log(lam))));
    fh = figure('position', [100 100 640 470]);
    plot(kk, gg, 'o-', 'color', [0 0.28 0.73], 'linewidth', 1.5, ...
         'markerfacecolor', [0 0.28 0.73], 'markersize', 4.5); hold on;
    plot(kk, g0, 's-', 'color', [0.78 0.12 0.12], 'linewidth', 1.5, ...
         'markerfacecolor', [0.78 0.12 0.12], 'markersize', 4.5);
    plot(kk([1 end]), gxg*[1 1], '--', 'color', [0 0.28 0.73], 'linewidth', 1);
    plot(kk([1 end]), gx0*[1 1], '--', 'color', [0.78 0.12 0.12], 'linewidth', 1);
    plot(kh, gh0, 'ko', 'markersize', 8, 'linewidth', 1.2);
    text(kh(1), gh0(1)-0.06, 'h = 1', 'horizontalalignment','center','fontsize',8);
    text(kh(2), gh0(2)+0.06, 'h = 2', 'horizontalalignment','center','fontsize',8);
    grid on; xlabel('predictor support size k (nearest even vertices)');
    ylabel('transform coding gain (dB)');
    ylim([min(g0)-0.14, max([gg gxg])+0.12]); xlim([0 kk(end)+1]);
    legend({'k-NN Wiener + Gouze', 'k-NN Wiener, U = 0', ...
            'exact + Gouze', 'exact, U = 0', 'hop-restricted (Table VII)'}, ...
            'location', 'southeast', 'fontsize', 8);
    legend boxoff;
    ax2 = axes('position', [0.30 0.245 0.30 0.27]);
    semilogy(kk, max(gx0 - g0, 1e-4), 's-', 'color', [0.78 0.12 0.12], ...
             'linewidth', 1.2, 'markerfacecolor', [0.78 0.12 0.12], ...
             'markersize', 3); grid on;
    set(ax2, 'fontsize', 7); ylim([5e-5 2]);
    set(ax2, 'ytick', [1e-4 1e-3 1e-2 1e-1 1]);
    xlabel('k', 'fontsize', 8); ylabel('G_{exact} - G(k)  (dB)', 'fontsize', 8);
    title('localization deficit', 'fontsize', 8, 'fontweight', 'normal');
    saveboth(fh, od, 'fig4_localization_sweep');
    fprintf('  fig4 k-sweep U=0: %s\n', mat2str(g0, 5));
    fprintf('  fig4 k-sweep Gz : %s\n', mat2str(gg, 5));
    fprintf('  fig4 hop pts: h=1 @ %.1f -> %.3f | h=2 @ %.1f -> %.3f (exact %.3f/%.3f, KLT %.3f)\n', ...
            kh(1), gh0(1), kh(2), gh0(2), gx0, gxg, gk);
catch err
    warning('Fig 4 skipped: %s', err.message);
end

% =====================================================================
% TABLE 7  (Laplacian learning vs sample size)
% =====================================================================
fprintf('\n== Table 7 (graph learning) ==\n');
G2 = sensor_graph_tests(30, 4);
delta = 0.5; Jt = G2.L + delta*eye(30);
Ct = chol(inv(Jt)).';
Ns = [1e3 5e3 1e4 4e4 1e5];
er = zeros(size(Ns));
for i = 1:numel(Ns)
    seed(100 + i);
    X  = Ct * randn(30, Ns(i));
    Sh = X*X.'/Ns(i);
    wl = learn_laplacian(Sh, G2.B, delta, 400);
    er(i) = norm(wl - G2.w) / norm(G2.w);
    fprintf('  S = %-7d rel err = %.3f\n', Ns(i), er(i));
end
f7 = fopen(fullfile(od, 'table7_learning.tex'), 'w');
fprintf(f7, '$\\|\\hat w - w\\|_2 / \\|w\\|_2$ & %.3f & %.3f & %.3f & %.3f & %.3f\\\\\n', er);
fclose(f7);

fprintf('\nAll outputs written to ./%s/\n', od);
d = dir(fullfile(od, '*'));
for i = 1:numel(d)
    if ~d(i).isdir, fprintf('  %s\n', d(i).name); end
end
end

% ======================================================================
% local helpers
% ======================================================================
function saveboth(fh, od, name)
print(fh, '-dpng', '-r300', fullfile(od, [name '.png']));
print(fh, '-depsc', fullfile(od, [name '.eps']));
end

function A = adjof(L)
W = -L; W(1:size(W,1)+1:end) = 0;
A = W > 1e-9;
end

function h0 = adap_orth(r, N, p)
if numel(r) < 4*N, r = [r(:); zeros(4*N - numel(r), 1)]; end
d = design_compaction_fb(r, N, p);
h0 = d.h0;
end

function s = fbw(d)
s = struct('h0', d.h0, 'h1', d.h1, 'g0', d.g0, 'g1', d.g1);
end

function seed(k)
try, rng(k); catch, rand('seed', k); randn('seed', k); end %#ok
end

function G = sensor_graph_demo(n, kseed)
% identical to demo_graph_benchmark (Tables 5-6, Figs 3-4)
seed(kseed);
for tryi = 1:50
    xy = rand(n, 2);
    Dm = sqrt((xy(:,1)-xy(:,1).').^2 + (xy(:,2)-xy(:,2).').^2);
    rad = 0.11 + 0.005*tryi;
    A = (Dm < rad) & ~eye(n);
    W = exp(-Dm.^2/(2*0.08^2)) .* A;
    L = diag(sum(W,2)) - W;
    ev = sort(eig(L));
    if ev(2) > 1e-8, break; end
end
G = struct('xy', xy, 'W', W, 'L', L, 'A', A);
end

function G = sensor_graph_tests(n, kseed)
% identical to run_tests3 (Table 7), includes edge incidence B
seed(kseed);
for tryi = 1:50
    xy = rand(n, 2);
    Dm = sqrt((xy(:,1)-xy(:,1).').^2 + (xy(:,2)-xy(:,2).').^2);
    rad = 0.18 + 0.01*tryi;
    A = (Dm < rad) & ~eye(n);
    W = exp(-Dm.^2 / (2*0.1^2)) .* A;
    L = diag(sum(W,2)) - W;
    ev = sort(eig(L));
    if ev(2) > 1e-8, break; end
end
[ii, jj] = find(triu(A, 1));
m = numel(ii);
B = zeros(n, m); ww = zeros(m, 1);
for k = 1:m
    B(ii(k), k) = 1; B(jj(k), k) = -1;
    ww(k) = W(ii(k), jj(k));
end
G = struct('xy', xy, 'W', W, 'L', L, 'B', B, 'w', ww, 'A', A);
end
