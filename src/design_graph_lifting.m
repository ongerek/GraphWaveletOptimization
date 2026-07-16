function out = design_graph_lifting(L, Sigma, e, opts)
%DESIGN_GRAPH_LIFTING One level of statistics-adapted lifting on a graph.
%
%   out = DESIGN_GRAPH_LIFTING(L, Sigma, e, opts)
%   L     : graph Laplacian (n x n);   Sigma : exact signal covariance;
%   e     : logical mask of 'even' (kept) vertices.
%
%     split:    x_e = x(e),  x_o = x(~e)
%     predict:  d = x_o - P x_e            (PR structural for ANY P, U)
%     update:   a = x_e + U d
%
%   Predict options (opts.predict):
%     'average'  p_ij = w_ij / sum_j w_ij over even neighbors (topology
%                only; the Narang-Ortega-style baseline)
%     'harmonic' P = -L_oo^{-1} L_oe  (harmonic extension; rows sum to 1
%                on a connected graph => d annihilates constants: one
%                graph vanishing moment, structurally)
%     'wiener'   MMSE predictor restricted to opts.hops-hop even
%                neighborhoods: per odd node i with support S,
%                p_i = Sigma(S,S)^{-1} Sigma(S,i)   (graph analog of the
%                np-tap constrained Wiener predictor of Module 2)
%     'exact'    P = Sigma_oe' ... full conditional mean
%                E[x_o|x_e] = Sigma_oe Sigma_ee^{-1} x_e  (dense; equals
%                -J_oo^{-1} J_oe for J = Sigma^{-1})
%
%   Update options (opts.update):
%     'none'     U = 0
%     'gouze'    minimize reconstruction MSE with the detail branch
%                discarded:  J(U) = E||Ud||^2 + E||d - PUd||^2, giving
%                (I + P'P) U = P' (unconstrained closed form, independent
%                of Cov(d)); with opts.uhops support restriction the
%                masked normal equations use Cov(d) exactly.
%
%   Output: P, U, e; a-branch covariance Sigma_a (EXACT:
%   a = (I-UP) x_e + U x_o); detail covariance Sigma_d; per-level
%   analysis/synthesis blocks T, S with [a; d] = T [x_e; x_o].

if nargin < 4, opts = struct(); end
ptype = getdef(opts, 'predict', 'wiener');
utype = getdef(opts, 'update', 'gouze');
hops  = getdef(opts, 'hops', 1);
uhops = getdef(opts, 'uhops', 0);       % 0 => unconstrained closed form

ie = find(e); io = find(~e);
ne = numel(ie); no = numel(io);
W  = -L; W(1:size(W,1)+1:end) = 0;

switch ptype
    case 'average'
        P = W(io, ie);
        rs = sum(P, 2);
        P  = P ./ max(rs, eps);
        P(rs == 0, :) = wiener_rows(Sigma, io(rs == 0), ie, ...
                                    hopmask(W, io(rs == 0), ie, 2));
    case 'harmonic'
        P = -L(io, io) \ L(io, ie);
    case 'exact'
        P = Sigma(io, ie) / Sigma(ie, ie);
    case 'wiener'
        P = wiener_rows(Sigma, io, ie, hopmask(W, io, ie, hops));
    otherwise
        error('unknown predict');
end

% exact second-order statistics of the detail branch
Ad = [-P, eye(no)];                     % d = Ad * [x_e; x_o]
Sxx = [Sigma(ie, ie), Sigma(ie, io); Sigma(io, ie), Sigma(io, io)];
Sd  = Ad * Sxx * Ad.';  Sd = (Sd + Sd.')/2;

switch utype
    case 'none'
        U = zeros(ne, no);
    case 'gouze'
        if uhops == 0
            U = (eye(ne) + P.'*P) \ P.';
        else
            U = gouze_masked(P, Sd, hopmask(W, ie, io, uhops));
        end
    otherwise
        error('unknown update');
end

T  = [eye(ne) - U*P, U; -P, eye(no)];   % [a; d] = T [x_e; x_o]
S  = [eye(ne), -U; P, eye(no) - P*U];   % exact inverse (lifting)
Aa = [eye(ne) - U*P, U];
Sa = Aa * Sxx * Aa.';  Sa = (Sa + Sa.')/2;

out = struct('P', P, 'U', U, 'e', e, 'T', T, 'S', S, ...
             'Sigma_a', Sa, 'Sigma_d', Sd, ...
             'sigma_d2', mean(diag(Sd)), 'ie', ie, 'io', io);
end

% ----------------------------------------------------------------------
function v = getdef(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end

function M = hopmask(W, rows, cols, h)
% M(i,j) true if cols(j) is within h hops of rows(i) in graph W
n = size(W, 1);
A = double(W > 0);
R = A;
for k = 2:h, R = double((R*A + R) > 0); end
M = R(rows, cols) > 0;
end

function P = wiener_rows(Sigma, io, ie, M)
% per-row support-constrained MMSE predictor
P = zeros(numel(io), numel(ie));
for i = 1:numel(io)
    S = find(M(i, :));
    if isempty(S), [~, S] = max(Sigma(ie, io(i)).^2); end
    cols = ie(S);
    P(i, S) = (Sigma(cols, cols) \ Sigma(cols, io(i))).';
end
end

function U = gouze_masked(P, Sd, M)
% minimize tr(U Sd U') + tr((I-PU) Sd (I-PU)') over supp(U) in M.
% vec-form Hessian H[(r,c),(r',c')] = (I+P'P)(r,r') * Sd(c,c'); only the
% supported entries are assembled (never the full Kronecker product).
ne = size(P, 2); no = size(P, 1);
Q  = eye(ne) + P.'*P;
idx = find(M(:));
[ri, ci] = ind2sub([ne, no], idx);
H  = Q(ri, ri) .* Sd(ci, ci);
bf = P.' * Sd;
u  = zeros(ne*no, 1);
u(idx) = H \ bf(idx);
U = reshape(u, ne, no);
end
