function [w, Lhat] = learn_laplacian(Shat, B, delta, niter)
%LEARN_LAPLACIAN Attractive-GMRF graph learning (simplified
%   Egilmez-Pavez-Ortega): given sample covariance Shat and candidate
%   edge incidence B (n x m, columns b_e = e_i - e_j), solve the CONVEX
%   problem
%       min_{w >= 0}  tr(J Shat) - log det J,
%       J = delta*I + sum_e w_e b_e b_e'   (Laplacian-constrained MLE)
%   by projected gradient with backtracking. Gradient:
%       df/dw_e = b_e' Shat b_e - b_e' J^{-1} b_e.
if nargin < 4, niter = 300; end
n = size(Shat, 1); m = size(B, 2);
w = ones(m, 1) * 0.1;
c1 = sum(B .* (Shat * B), 1).';         % b_e' Shat b_e (constant)
f  = obj(w, B, Shat, delta, n);
eta = 0.1;
for it = 1:niter
    J    = delta*eye(n) + B*diag(w)*B.';
    Ji   = inv(J);
    g    = c1 - sum(B .* (Ji * B), 1).';
    for bt = 1:30
        wn = max(0, w - eta*g);
        fn = obj(wn, B, Shat, delta, n);
        if fn <= f - 1e-12, break; end
        eta = eta/2;
    end
    if norm(wn - w) < 1e-10 * (1 + norm(w)), w = wn; break; end
    w = wn; f = fn; eta = eta * 1.5;
end
Lhat = B*diag(w)*B.';
end
function f = obj(w, B, Shat, delta, n)
J = delta*eye(n) + B*diag(w)*B.';
[R, p] = chol((J+J.')/2);
if p > 0, f = Inf; return; end
f = sum(sum(J .* Shat)) - 2*sum(log(diag(R)));
end
