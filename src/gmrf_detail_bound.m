function s2 = gmrf_detail_bound(r, n)
%GMRF_DETAIL_BOUND Minimum achievable detail variance: the conditional
%   variance of an interior odd sample given ALL even samples, by exact
%   Gaussian conditioning on an n x n Toeplitz section:
%       Cov(x_o | x_e) = (J_oo)^{-1},   J = R^{-1}.
%   This is the infinite-order predict-lifting bound (Markov blanket
%   collapses it to 2 taps for first-order Gauss-Markov inputs).
if nargin < 2, n = 513; end
R  = toeplitz(r(1:n));
J  = inv(R);
io = 2:2:n;
C  = inv(J(io, io));
s2 = C(round(end/2), round(end/2));
end
