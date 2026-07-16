function demo_graph_benchmark
%DEMO_GRAPH_BENCHMARK Statistics-adapted graph lifting: benchmarks.
%
%   Setting A (unification check): path graph, AR(1)-matched GMRF.
%     The graph pyramid's transform coding gain should match the
%     classical Module-2 Wiener-lifting coding gain (boundary effects
%     aside), and both sit below the GFT = KLT bound.
%
%   Setting B: random geometric sensor network (n = 200), attractive
%     GMRF x ~ N(0, (L + 0.1 I)^{-1}) -- a proxy for hydro-meteorological
%     fields. Predict variants x update variants, 3 levels, vs KLT.

% ---------------- A: path graph = classical case ----------------------
% path bootstrap: locate src/ relative to this file
if exist('design_compaction_fb', 'file') ~= 2
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
end

rho = 0.95; n = 256;
J = (diag([1, (1+rho^2)*ones(1,n-2), 1]) + ...
     diag(-rho*ones(1,n-1),1) + diag(-rho*ones(1,n-1),-1)) / (1-rho^2);
Sg = inv(J);
Lp = diag([1, 2*ones(1,n-2), 1]) - diag(ones(1,n-1),1) - diag(ones(1,n-1),-1);
fprintf('\n=== A: path graph, AR(1) rho=0.95 GMRF, 4 levels ===\n');
for pt = {'average', 'harmonic', 'wiener', 'exact'}
    pyr = graph_lifting_pyramid(Lp, Sg, 4, ...
          struct('predict', pt{1}, 'hops', 1, 'update', 'gouze'));
    fprintf('  graph lifting, predict=%-9s : %7.3f dB\n', ...
            pt{1}, tc_gain(pyr.T, pyr.S, Sg));
end
lam = eig(Sg);
fprintf('  GFT = KLT bound (AM/GM eigs)    : %7.3f dB\n', ...
        10*log10(mean(lam)/exp(mean(log(lam)))));
r  = acf_from_psd(@(w) psd_ar([1 -rho], 1-rho^2, w)); r = r/r(1);
Gc = coding_gain_bior_adapt(@(rr) fbw(design_wiener_lifting(rr, 2, 2)), r, 4);
fprintf('  Module-2 classical WL(2,2)      : %7.3f dB   (stationary, no boundary)\n', Gc);

% ---------------- B: sensor network -----------------------------------
n = 200;
G = sensor_graph_local(n, 11);
delta = 0.1;
Sg = inv(G.L + delta*eye(n));
Sg = Sg / mean(diag(Sg));
fprintf('\n=== B: sensor network n=%d (|E|=%d), GMRF delta=%.2f, 3 levels ===\n', ...
        n, nnz(triu(G.A,1)), delta);
fprintf('  %-34s %10s %10s\n', 'predict', 'U=none', 'U=gouze');
cfg = {'average', 0; 'harmonic', 0; 'wiener', 1; 'wiener', 2; 'exact', 0};
% fully localized variant: 1-hop predict AND 1-hop-masked update
pyl = graph_lifting_pyramid(G.L, Sg, 3, struct('predict','wiener', ...
      'hops',1,'update','gouze','uhops',1));

for i = 1:size(cfg, 1)
    nm = cfg{i,1}; h = cfg{i,2};
    lab = nm; if strcmp(nm,'wiener'), lab = sprintf('wiener %d-hop', h); end
    g = zeros(1,2); ut = {'none', 'gouze'};
    for k = 1:2
        pyr = graph_lifting_pyramid(G.L, Sg, 3, ...
              struct('predict', nm, 'hops', max(h,1), 'update', ut{k}));
        g(k) = tc_gain(pyr.T, pyr.S, Sg);
    end
    fprintf('  %-34s %10.3f %10.3f\n', lab, g);
end
fprintf('  %-34s %21.3f\n', 'wiener 1-hop + 1-hop update', tc_gain(pyl.T, pyl.S, Sg));
lam = eig(Sg);
fprintf('  %-34s %21.3f\n', 'GFT = KLT bound', ...
        10*log10(mean(lam)/exp(mean(log(lam)))));
pyr = graph_lifting_pyramid(G.L, Sg, 3, struct('predict','wiener','hops',1));
szs = cellfun(@(d) numel(d.io), pyr.lev);
fprintf('  detail counts per level: %s | coarse: %d\n', ...
        mat2str(szs), numel(pyr.coarse));
fprintf(['\nNotes: gains are Katto-Yasuda transform coding gains from the\n' ...
         'EXACT covariance. ''harmonic'' uses only the topology (P = -L_oo^{-1}L_oe,\n' ...
         'one graph vanishing moment); ''wiener h-hop'' is the support-constrained\n' ...
         'MMSE predictor (graph analog of Module 2 taps); ''exact'' is the full\n' ...
         'conditional mean. Coarsening: exact covariance propagation + Kron\n' ...
         'reduction of the Laplacian.\n']);
end

function s = fbw(d)
s = struct('h0', d.h0, 'h1', d.h1, 'g0', d.g0, 'g1', d.g1);
end

function G = sensor_graph_local(n, kseed)
try, rng(kseed); catch, rand('seed',kseed); randn('seed',kseed); end %#ok
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
