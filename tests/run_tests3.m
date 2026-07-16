function fails = run_tests3
%RUN_TESTS3 Verification harness for the graph lifting module.
% path bootstrap: locate src/ relative to this file
if exist('design_compaction_fb', 'file') ~= 2
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
end

tr = @(name, err, tol) fprintf('%-56s %.3e  [%s]\n', name, err, ok(err<tol));
fails = 0;

% T1: path graph + AR(1)-GMRF reproduces Module 2 exactly --------------
% Precision of AR(1): tridiagonal, J_ii ~ (1+rho^2), J_{i,i+1} ~ -rho.
rho = 0.95; n = 256;
J = diag([1, (1+rho^2)*ones(1,n-2), 1]) + ...
    diag(-rho*ones(1,n-1), 1) + diag(-rho*ones(1,n-1), -1);
J = J / (1 - rho^2);
Sg = inv(J);
Lp = diag([1, 2*ones(1,n-2), 1]) - diag(ones(1,n-1),1) - diag(ones(1,n-1),-1);
e  = graph_bipartition(Lp);
alt = false(n,1); alt(1:2:end) = true;
ebit = min(max(abs(double(e(:)) - double(alt))), ...
           max(abs(double(e(:)) - double(~alt))));
tr('T1 path bipartition = even/odd split', ebit, 0.5);
fails = fails + (ebit >= 0.5);
dl = design_graph_lifting(Lp, Sg, e, struct('predict','wiener','hops',1));
pin = dl.P(10, :); pin = pin(pin ~= 0);
err = max(abs(sort(pin) - rho/(1+rho^2)));
tr('T1 interior predictor taps = rho/(1+rho^2)', err, 1e-10);
fails = fails + (err >= 1e-10);
err = abs(dl.Sigma_d(10,10) - (1-rho^2)/(1+rho^2));
tr('T1 interior sigma_d^2 = (1-rho^2)/(1+rho^2)', err, 1e-10);
fails = fails + (err >= 1e-10);

% T2: Kron reduction of J = exact marginal precision -------------------
Jc = kron_reduce(J, e);                     % Schur complement of J
err = max(max(abs(inv(Jc) - Sg(e, e))));
tr('T2 Kron(J) = marginal precision (inv = Sigma_ee)', err, 1e-8);
fails = fails + (err >= 1e-8);
% Kron of a Laplacian is a Laplacian with nonneg weights
Lc = kron_reduce(Lp, e);
offmax = max(max(triu(Lc, 1)));
rsmax  = max(abs(sum(Lc, 2)));
tr('T2 Kron(L) Laplacian: offdiag<=0, rowsum=0', max(offmax, rsmax), 1e-10);
fails = fails + (max(offmax, rsmax) >= 1e-10);

% T3: pyramid PR + harmonic vanishing moment ---------------------------
G  = sensor_graph(120, 5);
Jg = G.L + 0.1*eye(120); Sgg = inv(Jg);
for pt = {'average', 'harmonic', 'wiener', 'exact'}
    pyr = graph_lifting_pyramid(G.L, Sgg, 3, struct('predict', pt{1}));
    err = max(max(abs(pyr.S * pyr.T - eye(120))));
    tr(sprintf('T3 pyramid PR, predict=%s', pt{1}), err, 1e-9);
    fails = fails + (err >= 1e-9);
end
dl = design_graph_lifting(G.L, Sgg, graph_bipartition(G.L), ...
                          struct('predict', 'harmonic'));
err = max(abs(dl.T(numel(dl.ie)+1:end, :) * ones(120,1)));
tr('T3 harmonic predict: d(const) = 0 (graph VM)', err, 1e-9);
fails = fails + (err >= 1e-9);

% T4: Gouze update optimality (perturbation) ---------------------------
dl = design_graph_lifting(G.L, Sgg, graph_bipartition(G.L), ...
        struct('predict', 'wiener', 'hops', 1, 'update', 'gouze'));
Jfun = @(U) trace(U*dl.Sigma_d*U.') + ...
    trace((eye(size(dl.P,1)) - dl.P*U)*dl.Sigma_d*(eye(size(dl.P,1)) - dl.P*U).');
J0 = Jfun(dl.U); worst = -Inf;
for k = 1:5
    E = zeros(size(dl.U)); E(randi(numel(E))) = 1e-4;
    worst = max(worst, max(J0 - Jfun(dl.U+E), J0 - Jfun(dl.U-E)));
end
tr('T4 Gouze update local optimality', max(worst, 0), 1e-10);
fails = fails + (worst >= 1e-10);

% T4b: masked Gouze with full support == unconstrained closed form ----
dlf = design_graph_lifting(G.L, Sgg, graph_bipartition(G.L), ...
        struct('predict','wiener','hops',1,'update','gouze','uhops',99));
dlu = design_graph_lifting(G.L, Sgg, graph_bipartition(G.L), ...
        struct('predict','wiener','hops',1,'update','gouze'));
err = max(max(abs(dlf.U - dlu.U)));
tr('T4b masked Gouze (full mask) = closed form', err, 1e-8);
fails = fails + (err >= 1e-8);

% T5: tc_gain of GFT/KLT = AM/GM of eigenvalues ------------------------
[V, D] = eig(Sgg); lam = diag(D);
gk = tc_gain(V.', V, Sgg);
gt = 10*log10(mean(lam)/exp(mean(log(lam))));
tr('T5 tc_gain(GFT) = AM/GM eigenvalue bound', abs(gk-gt), 1e-9);
fails = fails + (abs(gk-gt) >= 1e-9);

% T6: Laplacian learning recovers edge weights -------------------------
try, rng(3); catch, rand('seed',3); randn('seed',3); end %#ok
G2 = sensor_graph(30, 4);
n2 = 30; delta = 0.5;
Jt = G2.L + delta*eye(n2);
X  = chol(inv(Jt)).' * randn(n2, 40000);
Sh = X*X.'/40000;
[wl, ~] = learn_laplacian(Sh, G2.B, delta, 400);
err = norm(wl - G2.w) / norm(G2.w);
tr('T6 Laplacian learning: rel weight error (40k samp)', err, 0.12);
fails = fails + (err >= 0.12);

fprintf('\n%d test group failures.\n', fails);
end

function s = ok(b), if b, s='PASS'; else, s='FAIL'; end, end

function G = sensor_graph(n, kseed)
% connected random geometric graph with Gaussian kernel weights
try, rng(kseed); catch, rand('seed',kseed); randn('seed',kseed); end %#ok
for tryi = 1:50
    xy = rand(n, 2);
    Dm = sqrt((xy(:,1)-xy(:,1).').^2 + (xy(:,2)-xy(:,2).').^2);
    rad = 0.18 + 0.01*tryi;
    A = (Dm < rad) & ~eye(n);
    W = exp(-Dm.^2 / (2*0.1^2)) .* A;
    L = diag(sum(W,2)) - W;
    ev = sort(eig(L));
    if ev(2) > 1e-8, break; end          % connected
end
[ii, jj] = find(triu(A, 1));
m = numel(ii);
B = zeros(n, m); w = zeros(m,1);
for k = 1:m
    B(ii(k), k) = 1; B(jj(k), k) = -1;
    w(k) = W(ii(k), jj(k));
end
G = struct('xy', xy, 'W', W, 'L', L, 'B', B, 'w', w, 'A', A);
end
