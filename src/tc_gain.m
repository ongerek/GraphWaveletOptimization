function GdB = tc_gain(T, S, Sigma)
%TC_GAIN High-rate transform coding gain of a biorthogonal linear
%   transform pair y = T x, x = S y (Katto-Yasuda weighting by synthesis
%   column norms; reduces to AM/GM of eigenvalues for the KLT):
%       G = mean(diag(Sigma)) / geomean( var(y_i) * ||S(:,i)||^2 ).
vy = diag(T * Sigma * T.');
w  = sum(S.^2, 1).';
GdB = 10*log10( mean(diag(Sigma)) / exp(mean(log(vy .* w))) );
end
