function pyr = graph_lifting_pyramid(L, Sigma, Lev, opts)
%GRAPH_LIFTING_PYRAMID Multi-level statistics-adapted graph lifting.
%
%   pyr = GRAPH_LIFTING_PYRAMID(L, Sigma, Lev, opts) applies Lev levels of
%   DESIGN_GRAPH_LIFTING. At every level:
%     1. bipartition the CURRENT graph (greedy max-cut),
%     2. design predict/update from the CURRENT exact covariance,
%     3. coarsen: signal covariance propagated exactly through the
%        lifting step; graph topology by KRON REDUCTION of the current
%        Laplacian (with optional edge sparsification opts.kthresh) --
%        the Schur-complement coarsening that keeps the coarse signal a
%        GMRF on a legitimate graph.
%   Returns the full N x N analysis operator pyr.T (x -> [d_1; d_2; ...;
%   d_Lev; a_Lev] ordering: details first per level, coarse last), the
%   exact synthesis pyr.S (composed lifting inverses, NOT a matrix
%   inverse), per-level structures pyr.lev{j}, and pyr.perm bookkeeping.
if nargin < 4, opts = struct(); end
kthresh = 0; if isfield(opts, 'kthresh'), kthresh = opts.kthresh; end
n  = size(L, 1);
T  = eye(n); S = eye(n);
cur = 1:n;                              % global indices of current nodes
pyr.lev = {};
for j = 1:Lev
    nj = numel(cur);
    if nj < 4, break; end
    e  = graph_bipartition(L);
    dl = design_graph_lifting(L, Sigma, e, opts);
    % embed level transform into global operator
    Pi = [dl.ie(:); dl.io(:)];          % local permutation: evens first
    Tj = eye(n); Sj = eye(n);
    gl = cur(Pi);                       % global ids, permuted
    Tj(gl, gl) = dl.T;                  % rows/cols in permuted local order
    Sj(gl, gl) = dl.S;
    T = Tj * T;  S = S * Sj;
    % coarsen
    Sigma = dl.Sigma_a;
    L     = kron_reduce(L, e, kthresh);
    pyr.lev{j} = dl;
    cur = cur(dl.ie);
end
pyr.T = T; pyr.S = S; pyr.coarse = cur;
% reorder rows of T so output = [all details (fine->coarse); a_L]
ord = [];
kept = 1:n;
for j = 1:numel(pyr.lev)
    gl_o = kept(pyr.lev{j}.io);
    ord  = [ord, gl_o]; %#ok<AGROW>
    kept = kept(pyr.lev{j}.ie);
end
ord = [ord, kept];
pyr.T = pyr.T(ord, :);
pyr.S = pyr.S(:, ord);
pyr.order = ord;
end
