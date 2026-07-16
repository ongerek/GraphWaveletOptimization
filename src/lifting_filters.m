function [h0, h1, g0, g1, off] = lifting_filters(p, kp, u, ku)
%LIFTING_FILTERS Equivalent 2-channel filters of a predict/update
%   lifting scheme, extracted by probing (PR by construction):
%     d[n] = x[2n+1] - sum_i p_i x[2(n+kp_i)]
%     a[n] = x[2n]   + sum_j u_j d[n+ku_j]
%   Conventions:  a[m] = sum_k h0(k-off(1)+1) x[2m-k]   (same for h1),
%   synthesis  x[t] += g0(k-off(3)+1) a[m], k = t-2m    (same for g1).
n  = 128;  m0 = n/4;
h0 = zeros(1, 2*n); h1 = zeros(1, 2*n);
for t = 0:n-1
    x = zeros(n, 1); x(t+1) = 1;
    [a, d] = lift_fwd(x, p, kp, u, ku);
    k = 2*m0 - t;
    if abs(k) <= n-1
        h0(k+n+1) = a(m0+1);
        h1(k+n+1) = d(m0+1);
    end
end
av = zeros(n/2,1); dv = zeros(n/2,1);
av(m0+1) = 1;  g0 = readout(lift_inv(av, dv, p, kp, u, ku), n, m0);
av(:) = 0; dv(m0+1) = 1;
g1 = readout(lift_inv(av, dv, p, kp, u, ku), n, m0);
[h0, o0] = trim(h0, n); [h1, o1] = trim(h1, n);
[g0, o2] = trim(g0, n); [g1, o3] = trim(g1, n);
off = [o0 o1 o2 o3];
end

function [a, d] = lift_fwd(x, p, kp, u, ku)
s = x(1:2:end); t = x(2:2:end);
d = t;
for i = 1:numel(p), d = d - p(i)*circshift(s, -kp(i)); end
a = s;
for j = 1:numel(u), a = a + u(j)*circshift(d, -ku(j)); end
end

function x = lift_inv(a, d, p, kp, u, ku)
s = a;
for j = 1:numel(u), s = s - u(j)*circshift(d, -ku(j)); end
t = d;
for i = 1:numel(p), t = t + p(i)*circshift(s, -kp(i)); end
x = zeros(2*numel(s), 1);
x(1:2:end) = s; x(2:2:end) = t;
end

function g = readout(x, n, m0)
g = zeros(1, 2*n);
for t = 0:n-1
    k = t - 2*m0;
    if abs(k) <= n-1, g(k+n+1) = x(t+1); end
end
end

function [h, off] = trim(h, n)
i   = find(abs(h) > 1e-12);
off = i(1) - (n + 1);
h   = h(i(1):i(end));
end
