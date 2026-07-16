function demo_benchmark
%DEMO_BENCHMARK Signal-adapted orthonormal FB vs Daubechies and CDF 9/7.
%
%   Three WSS processes:
%     P1  AR(1), rho = 0.95                     (lowpass, "natural signal")
%     P2  AR(2), poles 0.975 e^{+-j 0.8 pi}     (high-frequency resonance)
%     P3  bimodal PSD: bumps at 0.15 pi / 0.70 pi + white floor
%   Compared on L = 4 level dyadic coding gain (dB):
%     db2, db4, db8, CDF 9/7, adapted N=16 p=2 (fixed filter),
%     adapted N=16 p=2 re-designed at every level (greedy),
%     and adapted N=16 p=0 (no vanishing-moment constraint).
%   Also reported: single-stage compaction sigma_a^2 vs the ideal
%   (infinite-order) compaction bound  (1/pi) int max(S(w), S(pi-w)) dw.

% path bootstrap: locate src/ relative to this file
if exist('design_compaction_fb', 'file') ~= 2
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
end

L = 4; N = 16;
K = 2^15; w = 2*pi*(0:K-1).'/K;

proc = struct('name', {}, 'S', {});
proc(1).name = 'AR(1) rho=0.95';
proc(1).S = @(w) psd_ar([1 -0.95], 1 - 0.95^2, w);
rho = 0.975; th = 0.8*pi; a2 = [1, -2*rho*cos(th), rho^2];
proc(2).name = 'AR(2) r=0.975 @0.8pi';
proc(2).S = @(w) psd_ar(a2, 1, w);
bump = @(w, w0, s) exp(-0.5*((w - w0)/s).^2) + exp(-0.5*((2*pi - w - w0)/s).^2);
proc(3).name = 'Bimodal + floor';
proc(3).S = @(w) bump(w, 0.15*pi, 0.06*pi) + 0.8*bump(w, 0.70*pi, 0.05*pi) + 0.01;

hdb = {db_filter(2), db_filter(4), db_filter(8)};
[b0, b1, bg0, bg1] = cdf97_filters();

fprintf('\n===== L = %d level dyadic coding gains (dB) =====\n\n', L);
fprintf('%-22s %7s %7s %7s %7s %8s %8s %8s\n', 'process', ...
        'db2', 'db4', 'db8', '9/7', 'adapt', 'adapt/L', 'adapt-p0');

for ip = 1:numel(proc)
    r = acf_from_psd(proc(ip).S, K);
    r = r / r(1);                                    % unit variance
    G = zeros(1, 7);
    for i = 1:3
        G(i) = coding_gain_ortho(@(rr) hdb{i}, r, L);
    end
    G(4) = coding_gain_bior(b0, b1, bg0, bg1, r, L);
    des  = design_compaction_fb(r, N, 2);
    G(5) = coding_gain_ortho(@(rr) des.h0, r, L);
    G(6) = coding_gain_ortho(@(rr) adapt(rr, N, 2), r, L);
    G(7) = coding_gain_ortho(@(rr) adapt(rr, N, 0), r, L);
    Sg = proc(ip).S(w);  Sg = Sg / mean(Sg);
    Gmax = -10*log10(exp(mean(log(Sg))));      % -10 log10 spectral flatness
    fprintf('%-22s %7.3f %7.3f %7.3f %7.3f %8.3f %8.3f %8.3f | bound %6.3f\n', ...
            proc(ip).name, G, Gmax);

    % single-stage compaction vs ideal bound
    h8 = hdb{3}; f8 = conv(h8, fliplr(h8)); c8 = numel(h8);
    sa_db = f8(c8)*r(1) + 2*(f8(c8+1:end) * r(2:c8));
    wf = linspace(0, pi, 2^14).';
    Sw = proc(ip).S(wf);  Sm = proc(ip).S(pi - wf);
    sa_id = trapz(wf, max(Sw, Sm)) / pi / (trapz(wf, Sw)/pi);  % normalized
    fprintf('%-22s single-stage sigma_a^2:  db8 %.4f | adapted %.4f | ideal %.4f\n\n', ...
            '', sa_db, des.sigma_a2, sa_id);
end

fprintf(['Notes: "adapt"   = N=16, p=2 VM, one LP design used at all levels;\n' ...
         '       "adapt/L" = re-optimized at each level (greedy stage-wise);\n' ...
         '       "adapt-p0"= no vanishing moments (pure compaction filter).\n' ...
         'For orthonormal FBs, sigma_a^2 in [0,2]; ideal bound is the\n' ...
         'infinite-order compaction limit (1/pi) int max(S(w),S(pi-w)) dw.\n\n']);

% ---- plot: bimodal case ----------------------------------------------
try
    r  = acf_from_psd(proc(3).S, K);  r = r / r(1);
    d2 = design_compaction_fb(r, N, 2);
    d0 = design_compaction_fb(r, N, 0);
    wf = linspace(0, pi, 1024).';
    Sw = proc(3).S(wf);  Sw = Sw / max(Sw) * 2;
    Fa2 = abs(freqz_local(d2.h0, wf)).^2;
    Fa0 = abs(freqz_local(d0.h0, wf)).^2;
    Fdb = abs(freqz_local(hdb{3}, wf)).^2;
    fig = figure('visible', 'off');
    plot(wf/pi, Sw, 'k-', wf/pi, Fdb, 'b--', ...
         wf/pi, Fa2, 'r-', wf/pi, Fa0, 'm-.', 'linewidth', 1.4);
    legend('PSD (scaled)', 'db8 |H_0|^2', 'adapted p=2 |H_0|^2', ...
           'adapted p=0 |H_0|^2', 'location', 'north');
    xlabel('\omega / \pi'); grid on;
    title('Bimodal PSD: adapted vs fixed lowpass responses');
    print(fig, '-dpng', '-r120', 'bimodal_responses.png');
    fprintf('Saved plot: bimodal_responses.png\n');
catch err
    fprintf('(plot skipped: %s)\n', err.message);
end
end

function h0 = adapt(r, N, p)
if numel(r) < 4*N, r = [r(:); zeros(4*N - numel(r), 1)]; end
d  = design_compaction_fb(r, N, p);
h0 = d.h0;
end

function H = freqz_local(h, w)
H = exp(-1i * w(:) * (0:numel(h)-1)) * h(:);
end
