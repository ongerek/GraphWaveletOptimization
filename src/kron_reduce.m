function Lc = kron_reduce(L, keep, thresh)
%KRON_REDUCE Kron reduction (Schur complement) of a graph Laplacian onto
%   the vertex subset 'keep' (logical or index):
%       Lc = L_kk - L_ke * L_ee^{-1} * L_ek
%   Theorem (Dorfler & Bullo, IEEE TCAS-I 2013): Lc is again a Laplacian
%   of a connected graph with NONNEGATIVE weights. If x ~ N(0, J^{-1})
%   with J = L + delta*I, the Kron reduction of J is EXACTLY the marginal
%   precision of x_keep (Schur complement = inverse of a principal
%   submatrix of the covariance) -- the coarse signal is again a GMRF on
%   a legitimate graph, which is what makes the multiscale recursion
%   statistically self-consistent.
%   thresh (optional): drop reduced edges with weight < thresh * max
%   weight (sparsification, cf. Shuman-Faraji-Vandergheynst pyramid);
%   the diagonal is recomputed so Lc stays an exact Laplacian.
if islogical(keep), keep = find(keep); end
n    = size(L, 1);
gone = setdiff(1:n, keep);
Lc   = L(keep, keep) - L(keep, gone) * (L(gone, gone) \ L(gone, keep));
Lc   = (Lc + Lc.')/2;
if nargin > 2 && thresh > 0
    W = -Lc; W(1:size(W,1)+1:end) = 0;
    W(W < thresh * max(W(:))) = 0;
    Lc = diag(sum(W, 2)) - W;
end
end
