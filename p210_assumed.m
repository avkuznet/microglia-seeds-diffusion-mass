function A = p210_assumed()
%P210_ASSUMED  How much the main results depend on each assumed input.
%
%   The inputs listed in the manuscript as assumed are neither measured nor
%   fitted.  Each is halved and doubled here, one at a time and without
%   recalibration, and three results are recomputed:
%
%     load4    control graft load at 4 weeks, % of area (calibration target)
%     irf8     Irf8-/- graft load relative to control at 4 weeks (the
%              out-of-sample prediction; measured 0.244 +/- 0.061)
%     share16  carrier share of the amyloid imported into the graft by
%              16 weeks, %
%
%   Because nothing is refitted, a change in load4 overstates the effect of
%   an input whose shift the calibration would absorb; the two ratios are
%   the more informative columns.
%
%   The function also reports how close the baseline run comes to the two
%   ceilings, nu_max and b_max, which decides whether they matter at all.
%
%   A = P210_ASSUMED() returns the table as a struct and prints it.

names = {'rho','nu_max','k_2','F_c','k_1','seedRatio','k_d','k_r','k_p','k_u', ...
         'b_max','F_h_pct','m_h','g_n','detRatio','D_s','lam_s','tau_s'};
label = {'rho','nu_max','k_2','m_c','k_1','sigma','k_d','k_r','k_p','k_u', ...
         'beta_max','F_h','m_h','g_n','mu_det/m_seed','D_s','lambda_s','tau_s'};
fac = [0.5 2];

[b4, bk, b16] = outputs('', 1);
A.names = label;  A.base = [b4 bk b16];
A.load4 = zeros(numel(names),2);  A.irf8 = A.load4;  A.share16 = A.load4;

fprintf('\n=== assumed inputs: halved and doubled, one at a time ===\n\n');
fprintf('  baseline: graft load at 4 wk %.4f %%, Irf8-/- ratio %.4f, carrier share at 16 wk %.4f %%\n\n', b4, bk, b16);
fprintf('  %-14s   load at 4 wk, %%     Irf8-/- ratio      carrier share, %%\n', 'input');
fprintf('  %-14s   x0.5      x2        x0.5     x2        x0.5     x2\n', '');
for k = 1:numel(names)
    for s = 1:2
        [A.load4(k,s), A.irf8(k,s), A.share16(k,s)] = outputs(names{k}, fac(s));
    end
    fprintf('  %-14s  %7.4f  %7.4f    %7.4f  %7.4f    %7.4f  %7.4f\n', label{k}, ...
            A.load4(k,1), A.load4(k,2), A.irf8(k,1), A.irf8(k,2), A.share16(k,1), A.share16(k,2));
end

% ---- how close the baseline comes to the two ceilings --------------------
%   The peaks fall between the few output times used above, so this needs a
%   run saved often.
sol = p210_solve(p210_params('baseline'), 0:0.25:112);
p = sol.p;
A.nuPeak   = max(sol.nu(:));
beta       = sol.b ./ max(sol.nu, 1e-12);
beta(sol.nu < 1e-3) = 0;                 % ignore cells with almost no carriers
A.betaPeak = max(beta(:)) / p.b_max;
fprintf('\n  largest microglial density in the baseline run: %.3f of homeostatic (ceiling nu_max = %.1f)\n', A.nuPeak, p.nu_max);
fprintf('  largest cargo per carrier in the baseline run:  %.4f of beta_max\n', A.betaPeak);
fprintf('  free level in the host: starts at m_h = %.2f, is %.4f after 1 day and %.4f after 4 weeks\n', ...
        p.m_h, sol.m(sol.t == 1, end), sol.m(sol.t == 28, end));

% ---- two cargo-handling rates set to values seen in cultured microglia ----
%   Majumdar et al. (2008): microglia degrade about 20% of internalized
%   fibrillar Abeta in 2-3 days and still hold over 80% of it after 3 days.
fprintf('\n  rates consistent with cultured microglia:\n');
A.culture = zeros(2,2);
cs = {'k_r', 0.05; 'k_d', 0.07};
for q = 1:2
    pb = p210_params('baseline');  pb.(cs{q,1}) = cs{q,2};
    pk = p210_params('irf8ko');    pk.(cs{q,1}) = cs{q,2};
    db = p210_post(p210_solve(pb, [0 28]));  dk = p210_post(p210_solve(pk, [0 28]));
    A.culture(q,:) = [db.loadPct(end), dk.loadPct(end)/db.loadPct(end)];
    fprintf('    %s = %.2f per day: graft load at 4 wk %.4f %%, Irf8-/- ratio %.4f\n', ...
            cs{q,1}, cs{q,2}, A.culture(q,1), A.culture(q,2));
end

% ---- an upper bound on k_p from the host ----------------------------------
%   Removing 95% of the carriers leaves the measured host load at 1.00 to
%   1.07 of control.  The ratio the model gives rises with k_p.
fprintf('\n  host load at 4 wk with 95%% of carriers removed, relative to control:\n');
kpv = [0.5 1 2 4];  A.hostRatio = zeros(size(kpv));
for q = 1:numel(kpv)
    p1 = p210_params('baseline');  p1.k_p = kpv(q) * p1.k_p;
    p2 = p1;  p2.nu0 = 0.05;
    d1 = p210_post(p210_solve(p1, [0 28]));  d2 = p210_post(p210_solve(p2, [0 28]));
    A.hostRatio(q) = d2.hostPct(end) / d1.hostPct(end);
    fprintf('    k_p x %.1f: %.4f\n', kpv(q), A.hostRatio(q));
end
end

% -------------------------------------------------------------------------
function [load4, irf8, share16] = outputs(name, f)
TS   = [0 1 28 112];
i4   = find(TS == 28);
solB = p210_solve(setp('baseline', name, f), TS);            dB = p210_post(solB);
dK   = p210_post(p210_solve(setp('irf8ko', name, f), [0 28]));
load4   = dB.loadPct(i4);
irf8    = dK.loadPct(end) / dB.loadPct(i4);
share16 = 100 * dB.carrShare(end);
end

% -------------------------------------------------------------------------
function p = setp(variant, name, f)
p = p210_params(variant);
if ~isempty(name)
    p.(name) = p.(name) * f;
    p = p210_derive(p);                  % F_c and F_h_pct feed derived values
end
end
