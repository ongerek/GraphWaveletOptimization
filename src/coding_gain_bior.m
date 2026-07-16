function [GdB, info] = coding_gain_bior(h0, h1, g0, g1, r, L)
%CODING_GAIN_BIOR L-level dyadic coding gain for a biorthogonal FB
%   (Katto-Yasuda weighting by equivalent synthesis basis norms):
%
%   G_L = sigma_x^2 / [ (sigma_aL^2 * wL)^{2^-L}
%                        * prod_j (sigma_dj^2 * w_j)^{2^-j} ],
%   w = ||equivalent synthesis filter||^2. Reduces to CODING_GAIN_ORTHO
%   when the FB is orthonormal (all w = 1).
fb = struct('h0', h0(:).', 'h1', h1(:).', 'g0', g0(:).', 'g1', g1(:).');
[GdB, info] = coding_gain_bior_adapt(@(rr) fb, r, L);
end
