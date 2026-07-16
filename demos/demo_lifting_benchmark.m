function demo_lifting_benchmark
%DEMO_LIFTING_BENCHMARK Wiener-lifting vs fixed wavelets vs Module 1.
%   Same three processes as DEMO_BENCHMARK, L = 4 levels.
%   WL(a,b)   = Wiener lifting, a predict taps, b update taps, fixed;
%   WL(a,b)/L = re-designed at every level.
%   Also: level-1 detail variance vs the GMRF conditional bound
%   (infinite-order predict-lifting limit).

% path bootstrap: locate src/ relative to this file
if exist('design_compaction_fb', 'file') ~= 2
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src'));
end

L = 4; K = 2^15;
proc(1).name = 'AR(1) rho=0.95';
proc(1).S = @(w) psd_ar([1 -0.95], 1-0.95^2, w);
rho = 0.975; th = 0.8*pi;
proc(2).name = 'AR(2) r=0.975 @0.8pi';
proc(2).S = @(w) psd_ar([1, -2*rho*cos(th), rho^2], 1, w);
bump = @(w,w0,s) exp(-0.5*((w-w0)/s).^2) + exp(-0.5*((2*pi-w-w0)/s).^2);
proc(3).name = 'Bimodal + floor';
proc(3).S = @(w) bump(w,0.15*pi,0.06*pi) + 0.8*bump(w,0.70*pi,0.05*pi) + 0.01;

[b0,b1,bg0,bg1] = cdf97_filters();
h8 = db_filter(8);
h8s = struct('h0',h8,'h1',fliplr(h8).*(-1).^(0:15), ...
             'g0',fliplr(h8),'g1',(-1).^(0:15).*h8);  % CQF synthesis
[l0,l1,lg0,lg1] = lifting_filters([.5;.5],0:1,[.25;.25],-1:0);
l53 = @(rr) struct('h0',l0,'h1',l1,'g0',lg0,'g1',lg1);

fprintf('\n===== L = %d level coding gains (dB), lifting module =====\n\n', L);
fprintf('%-22s %6s %6s %6s %8s %8s %8s %9s %9s\n', 'process', ...
        '5/3','9/7','db8','ad-p2/L','WL(2,2)','WL(6,6)','WL(6,6)/L','WLfree/L');
for ip = 1:3
    r = acf_from_psd(proc(ip).S, K); r = r/r(1);
    G = zeros(1,8);
    G(1) = coding_gain_bior_adapt(l53, r, L);
    G(2) = coding_gain_bior(b0,b1,bg0,bg1, r, L);
    G(3) = coding_gain_bior_adapt(@(rr) h8s, r, L);
    G(4) = coding_gain_ortho(@(rr) adap_orth(rr,16,2), r, L);
    w22  = design_wiener_lifting(r, 2, 2);
    G(5) = coding_gain_bior_adapt(@(rr) fbs(w22), r, L);
    w66  = design_wiener_lifting(r, 6, 6);
    G(6) = coding_gain_bior_adapt(@(rr) fbs(w66), r, L);
    G(7) = coding_gain_bior_adapt(@(rr) fbs(design_wiener_lifting(rr,6,6)), r, L);
    ofree = struct('dualvm', false, 'primalvm', false);
    G(8) = coding_gain_bior_adapt(@(rr) fbs(design_wiener_lifting(rr,20,6,ofree)), r, L);
    fprintf('%-22s %6.3f %6.3f %6.3f %8.3f %8.3f %8.3f %9.3f %9.3f\n', proc(ip).name, G);
    fprintf('%-22s level-1 sigma_d^2: WL(2,2) %.5f | WL(6,6) %.5f | GMRF bound %.5f\n\n', ...
            '', w22.sigma_d2, w66.sigma_d2, gmrf_detail_bound(r));
end
fprintf(['Notes: WL predict = constrained Wiener (sum p = 1, one dual VM);\n' ...
         'update = Gouze reconstruction-MSE optimum (sum u = 1/2, one primal\n' ...
         'VM). WLfree/L = 20-tap predict, 6-tap update, NO moment constraints,\n' ...
         're-designed per level. PR is structural for every variant; gains\n' ...
         'use Katto-Yasuda synthesis weighting.\n']);
end

function s = fbs(x)
if isstruct(x)
    s = struct('h0',x.h0,'h1',x.h1,'g0',x.g0,'g1',x.g1);
else
    error('bad');
end
end

function s = fbs4(h0,h1,g0,g1) %#ok<DEFNU>
s = struct('h0',h0,'h1',h1,'g0',g0,'g1',g1);
end

function h0 = adap_orth(r, N, p)
if numel(r) < 4*N, r = [r(:); zeros(4*N-numel(r),1)]; end
d = design_compaction_fb(r, N, p);
h0 = d.h0;
end
