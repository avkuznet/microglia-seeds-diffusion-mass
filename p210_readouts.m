function R = p210_readouts(which)
%P210_READOUTS  How far the results depend on the observation rule.
%
%   R = P210_READOUTS            all cases
%   R = P210_READOUTS([1 2])     selected cases
%
%   The graft load is a stained area, and the model has to say how a field
%   of deposits of mass F and number density c shows on a section.  The
%   paper uses Eq. (14), area = kap * c^(1/3) * F^(2/3) * vis(F/c).  This
%   routine changes that rule and nothing else, refits the same three
%   targets each time, and prints what the model then gives.
%
%   Cases
%     1  Eq. (14), as in the paper
%     2  section formula with the true thickness, T = 25 um:
%            area fraction = c (pi/4) d^2 (T + 2d/3) vis
%        Eq. (14) is its limit for deposits much thinner than the section.
%     3  planar section, T = 0: area proportional to mass, same visibility
%     4  Eq. (14) with a shallower visibility curve, exponent 2
%     5  Eq. (14) with a steeper visibility curve, exponent 6
%     6  Eq. (14), with intracellular cargo read through the same visibility
%        rule, one object per cell
%     7  Eq. (14), with ALL intracellular cargo read as stained area
%
%   Refitted in every case, as in p210_calibrate: k_a to the 80% depletion
%   of host free Abeta, perPercent to the 4-week graft load of 0.586%, and
%   m_seed to the BLZ945 load ratio of 0.256.  The carrier-field constants
%   (chi_s, fnu_graft) do not depend on the observation rule and are kept.
%   m_seed is searched only where the deposits of the reference host are
%   themselves at least half visible, detRatio * m_seed <= F_hRef / c_hRef;
%   beyond that the host would be declared at an area it does not show.
%   Where the BLZ945 target cannot be met inside that range the closest
%   attainable ratio is reported and marked.
%
%   Nothing here is used by the paper's main results; the stored calibration
%   p210_fitted.mat is read and not changed.

C = { ...
    {'Eq. (14), as in the paper',             struct()}, ...
    {'section formula, T = 25 um',            struct('readout','section','sectT',0.025)}, ...
    {'planar section, T = 0 (area ~ mass)',   struct('readout','section','sectT',0)}, ...
    {'Eq. (14), visibility exponent 2',       struct('visExp',2)}, ...
    {'Eq. (14), visibility exponent 6',       struct('visExp',6)}, ...
    {'Eq. (14) + cargo, same visibility',     struct('cargoRead',1)}, ...
    {'Eq. (14) + all cargo read',             struct('cargoRead',2)}};
if nargin < 1 || isempty(which), which = 1:numel(C); end

D = p210_data();
fprintf('\n=== observation rule: what the results depend on ===\n');
fprintf('  Each case is refitted to the same three targets: host free Abeta\n');
fprintf('  depleted 80%%, graft load %.3f%% at 4 weeks, BLZ945 load ratio %.3f.\n', ...
        D.irf8Load.val(1), D.deplete.rel);

R = cell(1, 0);
for q = which(:)'
    R{end+1} = oneCase(C{q}{1}, C{q}{2}, D); %#ok<AGROW>
end

fprintf('\n%-37s %8s %8s %8s %8s %8s %8s %8s %11s %9s\n', '', 'BLZ945', 'Irf8', 'aged', 'load', 'load', 'load', 'host', 'carrier', 'inert');
fprintf('%-37s %8s %8s %8s %8s %8s %8s %8s %11s %9s\n', '', 'ratio', 'ratio', 'ratio', '2 wk %', '4 wk %', '16 wk %', '4 wk %', 'share %', '16 wk %');
for q = 1:numel(R)
    r = R{q};
    if r.blzFitted, mark = ' '; else, mark = '*'; end
    fprintf('%-37s %7.3f%s %8.3f %8.3f %8.3f %8.3f %8.2f %8.2f %5.1f /%5.1f %9.1e\n', r.name, r.blz, mark, r.irf8, r.aged, ...
            r.load2, r.load4, r.load16, r.host4, r.share4, r.share16, r.inert16);
end
fprintf('%-37s %8.3f %8.3f %8.3f %8.3f %8s %8.2f\n', 'measured', D.deplete.rel, D.irf8Load.rel, D.aged.rel, ...
        D.graftTime.val(1), '0.586', D.graftTime.val(3));
fprintf('  measured s.e.m.: BLZ945 %.3f, Irf8 %.3f, aged %.3f; loads at 2 and 16 weeks are from\n', ...
        D.deplete.relErr, D.irf8Load.relErr, D.aged.relErr);
fprintf('  another cohort (ED Fig. 1c), %.3f +/- %.3f and %.2f +/- %.2f.\n', ...
        D.graftTime.val(1), D.graftTime.err(1), D.graftTime.val(3), D.graftTime.err(3));
fprintf('  * target not attainable; closest ratio shown.  Carrier share is at 4 / 16 weeks.\n');

fprintf('\n%-37s %10s %10s %10s %14s %12s\n', '', 'k_a', 'perPercent', 'm_seed', 'onset rim/ctr', 'Irf8 rim');
for q = 1:numel(R)
    r = R{q};
    fprintf('%-37s %10.4f %10.2f %10.3g %6.2f /%5.2f %12.2f\n', r.name, r.k_a, r.pp, r.m_seed, r.onset(1), r.onset(3), r.onsetK);
end
fprintf('  onsets in weeks, of the deposits (cargo is not included in the onset).\n');

for q = 1:numel(R)
    r = R{q};
    if ~isempty(r.diam)
        fprintf('\n  %s: deposit diameters that the section geometry implies, um\n', r.name);
        fprintf('    reference host at transplantation %.1f;  half-visibility %.1f\n', r.diam.ref, r.diam.det);
        fprintf('    host   at 2, 4, 16 weeks   %5.1f %5.1f %5.1f\n', r.diam.host);
        fprintf('    graft  at 2, 4, 16 weeks   %5.1f %5.1f %5.1f\n', r.diam.graft);
        fprintf('    Eq. (14) assumes deposits much thinner than the %.0f um section.\n', 1e3 * 0.025);
    end
end
end

% =========================================================================
function r = oneCase(name, ov, D)
NC = 400;  TS = [0 7 14 21 28];  T16 = [0 14 28 112];
L0 = load('p210_fitted.mat');
ka = L0.k_a;  pp = L0.perPercent;  ms = L0.m_seed;
TGT_PCT = D.irf8Load.val(1);  TGT_BLZ = D.deplete.rel;
fitted = true;
t0 = tic;
if ~isempty(fieldnames(ov))
    for outer = 1:16
        fa = @(k) depl(mk('baseline',ov,k,pp,ms,NC), TS) - 0.80;
        ka = bis(fa, 1e-3, 5, 1e-4);
        ppOld = pp;
        for it = 1:40
            L  = loadOf(mk('baseline',ov,ka,pp,ms,NC), TS);
            np = pp * (L / TGT_PCT)^0.75;  rel = abs(np - pp) / np;  pp = np;
            if rel < 5e-5, break; end
        end
        p0   = mk('baseline',ov,ka,pp,ms,NC);
        msHi = p0.F_hRef / (p0.c_hRef * p0.detRatio);     % reference host half visible
        fm = @(m) loadOf(mk('depleted',ov,ka,pp,m,NC), TS) / max(loadOf(mk('baseline',ov,ka,pp,m,NC), TS), eps) - TGT_BLZ;
        [msNew, fitted] = bis(fm, 1e-5, msHi, 1e-4);
        relm = abs(msNew - ms) / msNew;  ms = msNew;
        if abs(pp - ppOld)/pp < 2e-4 && relm < 5e-4, break; end
    end
    if ~fitted
        % m_seed sits at the end of its range: bring k_a and perPercent to it
        for outer = 1:6
            fa = @(k) depl(mk('baseline',ov,k,pp,ms,NC), TS) - 0.80;
            ka = bis(fa, 1e-3, 5, 1e-4);
            ppOld = pp;
            for it = 1:40
                L  = loadOf(mk('baseline',ov,ka,pp,ms,NC), TS);
                np = pp * (L / TGT_PCT)^0.75;  rel = abs(np - pp) / np;  pp = np;
                if rel < 5e-5, break; end
            end
            if abs(pp - ppOld)/pp < 2e-4, break; end
        end
    end
end
pB = mk('baseline',ov,ka,pp,ms,NC);
sB = p210_solve(pB, T16);  dB = p210_post(sB);
dK = p210_post(p210_solve(mk('irf8ko',   ov,ka,pp,ms,NC), T16));
dD = p210_post(p210_solve(mk('depleted', ov,ka,pp,ms,NC), T16));
dA = p210_post(p210_solve(mk('aged',     ov,ka,pp,ms,NC), T16));
dN = p210_post(p210_solve(mk('nocarrier',ov,ka,pp,ms,NC), T16));
r.name = name;  r.k_a = ka;  r.pp = pp;  r.m_seed = ms;  r.blzFitted = fitted;
r.load2 = dB.loadPct(2);  r.load4 = dB.loadPct(3);  r.load16 = dB.loadPct(4);
r.host4 = dB.hostPct(3);
r.irf8  = dK.loadPct(3) / dB.loadPct(3);
r.blz   = dD.loadPct(3) / dB.loadPct(3);
r.aged  = dA.loadPct(3) / dB.loadPct(3);
r.share4  = 100 * dB.carrShare(3);  r.share16 = 100 * dB.carrShare(4);
r.inert16 = dN.loadPct(4);
ok = ~isnan(dB.onsetT);
if any(ok), r.onset = [min(dB.onsetT(ok)) median(dB.onsetT(ok)) max(dB.onsetT(ok))] / 7; else, r.onset = [NaN NaN NaN]; end
okK = ~isnan(dK.onsetT);
if any(okK), r.onsetK = min(dK.onsetT(okK)) / 7; else, r.onsetK = NaN; end
r.diam = [];
if isfield(pB,'d_ref')
    G = sB.G;
    r.diam.ref = 1e3 * pB.d_ref;
    r.diam.det = 1e3 * pB.d_ref * (pB.detRatio * pB.m_seed / pB.mu_ref)^(1/3);
    for k = 2:4
        F = sB.F(k,:)';  c = max(sB.c(k,:)', pB.c_min);
        mg = sum(G.Vc(G.inGraft) .* F(G.inGraft)) / sum(G.Vc(G.inGraft) .* c(G.inGraft));
        mh = sum(G.Vc(G.inHost)  .* F(G.inHost))  / sum(G.Vc(G.inHost)  .* c(G.inHost));
        r.diam.graft(k-1) = 1e3 * pB.d_ref * (mg / pB.mu_ref)^(1/3);
        r.diam.host(k-1)  = 1e3 * pB.d_ref * (mh / pB.mu_ref)^(1/3);
    end
end
fprintf('  %-37s done in %3.0f s\n', name, toc(t0));
end

function p = mk(variant, ov, ka, pp, ms, NC)
p = p210_params(variant);
p.N = NC;  p.k_a = ka;  p.perPercent = pp;  p.m_seed = ms;
fn = fieldnames(ov);
for q = 1:numel(fn), p.(fn{q}) = ov.(fn{q}); end
p = p210_derive(p);
end

function L = loadOf(p, TS)
d = p210_post(p210_solve(p, TS));
L = d.loadPct(end);
end

function dep = depl(p, TS)
sol = p210_solve(p, TS);  G = sol.G;
dep = 1 - sum(G.Vc(G.inHost) .* sol.m(end,G.inHost)') / G.Vhost;
end

function [x, ok] = bis(f, lo, hi, tol)
flo = f(lo);  fhi = f(hi);  ok = true;
if abs(flo) < tol, x = lo; return; end
if abs(fhi) < tol, x = hi; return; end
if flo * fhi > 0
    % no sign change at the ends: look for one inside, else take the closest point
    xs = logspace(log10(lo), log10(hi), 9);  fs = arrayfun(f, xs);
    k  = find(fs(1:end-1) .* fs(2:end) <= 0, 1);
    if isempty(k), ok = false;  [~, j] = min(abs(fs));  x = xs(j);  return; end
    lo = xs(k);  hi = xs(k+1);  flo = fs(k);
end
for it = 1:45
    x  = sqrt(lo * hi);  fx = f(x);
    if abs(fx) < tol, break; end
    if fx * flo > 0, lo = x; flo = fx; else, hi = x; end
end
end
