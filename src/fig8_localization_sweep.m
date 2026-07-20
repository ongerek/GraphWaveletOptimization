function fig8_localization_sweep
%FIG8_LOCALIZATION_SWEEP Manuscript Fig. 8 (redesigned): localization /
%   performance trade-off of the support-restricted Wiener predictor on
%   the sensor network of Table VII.
%
%   Sweeps the per-vertex support size k (the k nearest even vertices in
%   weighted shortest-path distance; predict = 'wiener_knn' in
%   DESIGN_GRAPH_LIFTING), with and without the Gouze update, against
%   the exact-conditional-mean asymptotes. The hop-restricted designs of
%   Table VII are overlaid as open circles at their mean level-1 support
%   sizes. Inset: localization deficit G_exact - G(k) (update off) on a
%   log scale -- the graph-domain mirror of the tap-support sweep of
%   Table V.
%
%   Uses the same sensor graph as demo_graph_benchmark / runme
%   (n = 200, seed 11). Outputs fig8_localization_sweep.{png,eps,pdf}
%   to ./paper_out/ and echoes every number to the console.

od = 'paper_out';
if ~exist(od, 'dir'), mkdir(od); end

% ---- graph and GMRF (identical to runme / demo_graph_benchmark) ------
n  = 200;
Gs = sensor_graph_demo(n, 11);
Ss = inv(Gs.L + 0.1*eye(n));  Ss = Ss / mean(diag(Ss));
lam = eig(Ss);
gk  = 10*log10(mean(lam) / exp(mean(log(lam))));
fprintf('sensor graph: |E| = %d, KLT bound %.3f dB\n', ...
        nnz(triu(Gs.A, 1)), gk);

% ---- exact-conditional-mean asymptotes -------------------------------
px0 = graph_lifting_pyramid(Gs.L, Ss, 3, struct('predict','exact','update','none'));
pxg = graph_lifting_pyramid(Gs.L, Ss, 3, struct('predict','exact','update','gouze'));
gx0 = tc_gain(px0.T, px0.S, Ss);
gxg = tc_gain(pxg.T, pxg.S, Ss);

% ---- k-NN support sweep ----------------------------------------------
kk = [1 2 3 4 6 8 10 14 18 24];
g0 = nan(size(kk));  gg = nan(size(kk));
for i = 1:numel(kk)
    p0 = graph_lifting_pyramid(Gs.L, Ss, 3, struct('predict','wiener_knn', ...
         'knn', kk(i), 'update', 'none'));
    pg = graph_lifting_pyramid(Gs.L, Ss, 3, struct('predict','wiener_knn', ...
         'knn', kk(i), 'update', 'gouze'));
    g0(i) = tc_gain(p0.T, p0.S, Ss);
    gg(i) = tc_gain(pg.T, pg.S, Ss);
end

% ---- hop-restricted operating points (U = 0), mean level-1 support ---
e1 = graph_bipartition(Gs.L);
A  = double(Gs.W > 0);  A(1:n+1:end) = 0;
kh = zeros(1, 2);  gh0 = zeros(1, 2);
R = A;
for h = 1:2
    if h == 2, R = double((R*A + R) > 0); end
    kh(h) = mean(sum(R(~e1, e1) > 0, 2));
    ph = graph_lifting_pyramid(Gs.L, Ss, 3, struct('predict','wiener', ...
         'hops', h, 'update', 'none'));
    gh0(h) = tc_gain(ph.T, ph.S, Ss);
end

% ---- console dump (fill the caption from these) ----------------------
fprintf('k-sweep U=0 : %s\n', mat2str(g0, 5));
fprintf('k-sweep Gz  : %s\n', mat2str(gg, 5));
fprintf('asymptotes  : exact U=0 %.3f | exact Gouze %.3f\n', gx0, gxg);
fprintf('hop points  : h=1 at kbar=%.1f -> %.3f | h=2 at kbar=%.1f -> %.3f\n', ...
        kh(1), gh0(1), kh(2), gh0(2));
fprintf('deficits U=0: %s\n', mat2str(gx0 - g0, 3));
fprintf('PR check    : %.1e\n', norm(px0.S*px0.T - eye(n), 'fro'));

% ---- figure ----------------------------------------------------------
fh = figure('visible', 'off', 'position', [100 100 640 470]);
plot(kk, gg, 'o-', 'color', [0 0.28 0.73], 'linewidth', 1.5, ...
     'markerfacecolor', [0 0.28 0.73], 'markersize', 4.5); hold on
plot(kk, g0, 's-', 'color', [0.78 0.12 0.12], 'linewidth', 1.5, ...
     'markerfacecolor', [0.78 0.12 0.12], 'markersize', 4.5);
plot(kk([1 end]), gxg*[1 1], '--', 'color', [0 0.28 0.73], 'linewidth', 1);
plot(kk([1 end]), gx0*[1 1], '--', 'color', [0.78 0.12 0.12], 'linewidth', 1);
plot(kh, gh0, 'ko', 'markersize', 8, 'linewidth', 1.2);
text(kh(1), gh0(1)-0.06, 'h = 1', 'horizontalalignment','center','fontsize',8);
text(kh(2), gh0(2)+0.06, 'h = 2', 'horizontalalignment','center','fontsize',8);
grid on
xlabel('predictor support size k (nearest even vertices)');
ylabel('transform coding gain (dB)');
ylim([min(g0)-0.14, max([gg gxg])+0.12]);  xlim([0 kk(end)+1]);
legend({'k-NN Wiener + Gouze', 'k-NN Wiener, U = 0', ...
        'exact + Gouze', 'exact, U = 0', 'hop-restricted (Table VII)'}, ...
        'location', 'southeast', 'fontsize', 8);
legend boxoff
ax2 = axes('position', [0.30 0.245 0.30 0.27]);
semilogy(kk, max(gx0 - g0, 1e-4), 's-', 'color', [0.78 0.12 0.12], ...
         'linewidth', 1.2, 'markerfacecolor', [0.78 0.12 0.12], ...
         'markersize', 3); grid on
set(ax2, 'fontsize', 7);  ylim([5e-5 2]);
set(ax2, 'ytick', [1e-4 1e-3 1e-2 1e-1 1]);
xlabel('k', 'fontsize', 8);
ylabel('G_{exact} - G(k)  (dB)', 'fontsize', 8);
title('localization deficit', 'fontsize', 8, 'fontweight', 'normal');

print(fh, '-dpng', '-r300', fullfile(od, 'fig8_localization_sweep.png'));
print(fh, '-depsc',         fullfile(od, 'fig8_localization_sweep.eps'));
print(fh, '-dpdf',          fullfile(od, 'fig8_localization_sweep.pdf'));
fprintf('figure written to %s/\n', od);
end

% ----------------------------------------------------------------------
function seed(k)
try, rng(k); catch, rand('seed', k); randn('seed', k); end %#ok
end

function G = sensor_graph_demo(n, kseed)
% identical to runme / demo_graph_benchmark (Tables VI-VII graph)
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
