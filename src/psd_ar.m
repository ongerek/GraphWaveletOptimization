function S = psd_ar(a, s2e, w)
%PSD_AR PSD of AR process  x[n] = -a(2)x[n-1] - ... + e[n]:
%   S(w) = s2e / |A(e^{jw})|^2,  a = [1 a1 a2 ...].
A = exp(-1i * w(:) * (0:numel(a)-1)) * a(:);
S = s2e ./ abs(A).^2;
end
