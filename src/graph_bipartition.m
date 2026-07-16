function e = graph_bipartition(L)
%GRAPH_BIPARTITION Even/odd vertex split for graph lifting: greedy
%   weighted MAX-CUT with spectral initialization (sign pattern of the
%   eigenvector of the largest Laplacian eigenvalue = the "highest graph
%   frequency", which alternates like (-1)^n on a path and recovers the
%   classical even/odd split there). Returns logical mask e (even set).
n = size(L, 1);
W = -L; W(1:n+1:end) = 0;              % weight matrix
[V, D] = eig((L + L.')/2);
[~, i] = max(diag(D));
s = sign(V(:, i)); s(s == 0) = 1;
% greedy 1-swap sweeps: move node to the side maximizing cut weight
for pass = 1:50
    changed = false;
    for v = 1:n
        gain = s(v) * (W(v, :) * s);   % >0 => flipping v increases cut
        if gain > 1e-12
            s(v) = -s(v); changed = true;
        end
    end
    if ~changed, break; end
end
e = (s > 0);
if all(e) || ~any(e), e(1:2:n) = true; e(2:2:n) = false; end
% convention: 'even' = larger side (kept/coarse set)
if nnz(e) < n/2, e = ~e; end
end
