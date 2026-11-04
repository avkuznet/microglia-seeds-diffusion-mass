function out = p210_calibrate(verbose)
%P210_CALIBRATE  Fit the free constants to the published Source Data.
%
%   Target                                                       Parameter
%   ---------------------------------------------------------   ----------
%   graft microglial density 751 /mm^2 at 2 weeks                 chi_s
%   graft microglial density 517 /mm^2 at 4 weeks                 fnu_graft
%     (Source Data Fig. 1d).  The carrier field is an input to the amyloid
%     problem, not a prediction of it, and it has two published values: the
%     acute injury signal drives the density above the graft's capacity and
%     it then relaxes back to it.  Fitting only one forces a monotone rise.
%   host interstitial Abeta depleted 80% at 4 weeks               k_a
%     (Cirrito et al. 2003; Hong et al. 2011)
%   graft amyloid load 0.586% of area at 4 weeks                  perPercent
%     (Source Data Fig. 2h, Irf8+/- cohort).  A unit conversion only.
%   BLZ945 graft load is 0.256 of vehicle                         m_seed
%     (Source Data ED Fig. 7f,g).  The deposit mass one carrier release
%     leaves behind sets how many deposits a given carrier flux starts, and
%     with it how steeply the load follows the carrier number.  Fitted to
%     the depletion arm, which changes the carrier number and nothing else.
%     Uptake per aged cell is NOT fitted: it is the ex vivo proxy 0.494 of
%     ED Fig. 4d.  Nothing in the aged arm is fitted at all; see below.
%   microglia at a laser lesion rise 1.82-fold in 24 h            Q0les
%     (Source Data Fig. 3e)
%
%   SIX constants, not seven.  The aged arm no longer fits a release rate:
%   its uptake per cell is the ex vivo proxy 0.494 of ED Fig. 4d and its
%   host load is declared, so its graft load is a prediction too,
%   conditional on that proxy.
%
%   The Irf8-/- arm is fitted to NOTHING.  Its graft carrier density is
%   imposed at the measured 0.133 of control (Fig. 2g) and nothing else is
%   changed: motility is NOT scaled, because Fig. 2d,e measures the travel
%   of microglial PROCESSES and the model's drift is translocation of the
%   SOMA.  The amyloid carried per cell (ED Fig. 6h, 0.663 of control) is a
%   standing content rather than an uptake rate and is not charged to the
%   model by default; p.useIrf8Cargo explores charging it.  The graft
%   amyloid load is therefore a prediction, as is the host deposit density,
%   which ED Fig. 6b reports unchanged.
%
%   Results are written to p210_fitted.mat and read back by p210_params.

if nargin < 1, verbose = true; end

D = p210_data();
TGT_DEPL = 0.80;
TGT_PCT  = D.irf8Load.val(1);      % 0.586 % area, Source Data Fig. 2h
TGT_BLZ  = D.deplete.rel;          % 0.256
TGT_NU   = D.lesNu.rel24;          % 1.82
TS   = [0 7 14 21 28];
NC   = 240;

perPercent = 34.49;  k_a = 0.282;  m_seed = 4.62e-4;
chi_s = 8.354e-2;    fnu_graft = 0.2739;  Q0les = 107.5;

% ---- 1. the carrier field, from its own two measurements ----------------
if verbose, fprintf('\n--- carrier field: chemotaxis and graft capacity ---\n'); end
for it = 1:6
    fcs   = @(cs) nuAt(14, NC, TS, k_a, perPercent, m_seed, cs, fnu_graft) - D.graftNu.val(1);
    chi_s = bisect(fcs, 1e-3, 1.0, 5, 20, verbose, 'chi_s');
    ffn   = @(fn) nuAt(28, NC, TS, k_a, perPercent, m_seed, chi_s, fn) - D.graftNu.val(2);
    fnuNew = bisect(ffn, 0.05, 1.5, 5, 20, verbose, 'fnu_g');
    rel = abs(fnuNew - fnu_graft) / fnuNew;
    fnu_graft = fnuNew;
    if verbose
        fprintf('  it %d: chi_s = %.4g  fnu_graft = %.4g  (change %.2g)\n', it, chi_s, fnu_graft, rel);
    end
    if rel < 5e-3, break; end
end

% ---- 2-4.  k_a, perPercent and m_seed -----------------------------------
%   All three are coupled.  perPercent sets the model units of the declared
%   host deposit area and therefore how fast the host consumes monomer, to
%   which k_a is fitted; m_seed sets the deposit number in the graft and
%   therefore the graft load, to which perPercent is fitted.  They are
%   brought to consistency together.
for outer = 1:8
    if verbose, fprintf('\n--- outer %d: host monomer depletion ---\n', outer); end
    fa  = @(ka) depletion(NC, TS, ka, perPercent, m_seed, chi_s, fnu_graft) - TGT_DEPL;
    k_a = bisect(fa, 1e-4, 50, 2e-3, 24, verbose, 'k_a');

    if verbose, fprintf('--- outer %d: histological unit conversion ---\n', outer); end
    ppOld = perPercent;
    for it = 1:10
        L = loadOf('baseline', NC, TS, k_a, perPercent, m_seed, chi_s, fnu_graft);
        if L <= 0, break; end
        % the readout goes as kap*F^(2/3) and kap ~ perPercent^(-2/3) at
        % fixed deposit mass, so a step of (L/target)^(3/2) is close to
        % Newton; damp it, because the graft mass moves with F_h as well
        newPP = perPercent * (L / TGT_PCT)^(0.75);
        rel   = abs(newPP - perPercent) / newPP;
        perPercent = newPP;
        if verbose
            fprintf('  it %2d: perPercent = %.5g   graft load %.4f%%  (change %.2g)\n', ...
                    it, perPercent, L, rel);
        end
        if rel < 1e-3, break; end
    end

    if verbose, fprintf('--- outer %d: deposit mass per release event ---\n', outer); end
    fm = @(ms) blzRatio(NC, TS, k_a, perPercent, ms, chi_s, fnu_graft) - TGT_BLZ;
    msNew = bisect(fm, 1e-4, 3e0, 3e-3, 24, verbose, 'm_seed');
    relm  = abs(msNew - m_seed) / msNew;
    m_seed = msNew;

    if verbose
        fprintf('  outer %d: perPercent %.5g (change %.2g), m_seed %.5g (change %.2g)\n', ...
                outer, perPercent, abs(perPercent-ppOld)/perPercent, m_seed, relm);
    end
    if abs(perPercent - ppOld)/perPercent < 3e-3 && relm < 5e-3, break; end
end


% ---- 5. acute lesion signal strength ------------------------------------
if verbose, fprintf('\n--- acute lesion signal ---\n'); end
fq    = @(q) lesionRatio(q, perPercent, chi_s, m_seed) - TGT_NU;
Q0les = bisect(fq, 10, 2e4, 2e-2, 24, verbose, 'Q0les');

f_chi = 0.0;   % retired: the Irf8 arm no longer fits anything
save('p210_fitted.mat', 'k_a', 'perPercent', 'm_seed', 'Q0les', ...
     'chi_s', 'fnu_graft', 'f_chi');

% ---- report --------------------------------------------------------------
A = struct('NC',NC,'TS',TS,'k_a',k_a,'pp',perPercent,'ms',m_seed, ...
           'cs',chi_s,'fn',fnu_graft);
[Lw, Hw, Nw] = fullOf('baseline',  A);
[Lk, ~ , Nk] = fullOf('irf8ko',    A);
[Ld, ~ , ~ ] = fullOf('depleted',  A);

out.k_a        = k_a;
out.perPercent = perPercent;
out.m_seed     = m_seed;
out.chi_s      = chi_s;
out.fnu_graft  = fnu_graft;
out.Q0les      = Q0les;
out.depl       = depletion(NC, TS, k_a, perPercent, m_seed, chi_s, fnu_graft);
out.load_wt    = Lw;
out.hostPct    = Hw;
out.blz        = Ld / max(Lw, eps);
out.agedRatio  = agedRatio(NC, TS, k_a, perPercent, m_seed, chi_s, fnu_graft);
out.nuRatio    = lesionRatio(Q0les, perPercent, chi_s, m_seed);
out.ratio      = Lk / max(Lw, eps);          % PREDICTION
out.nuRatio_ko = Nk / max(Nw, eps);          % PREDICTION
out.nu2        = nuAt(14, NC, TS, k_a, perPercent, m_seed, chi_s, fnu_graft);
out.nu4        = nuAt(28, NC, TS, k_a, perPercent, m_seed, chi_s, fnu_graft);
out.gamma      = log(out.blz) / log(D.deplete.nuRel);

if verbose
    fprintf('\n=== fitted constants ===\n');
    fprintf('  chi_s      = %.5g         graft nu at 2 wk %.0f (target %.0f)\n', chi_s, out.nu2, D.graftNu.val(1));
    fprintf('  fnu_graft  = %.5g         graft nu at 4 wk %.0f (target %.0f)\n', fnu_graft, out.nu4, D.graftNu.val(2));
    fprintf('  k_a        = %.5g 1/day   host depletion %.3f (target %.2f)\n', k_a, out.depl, TGT_DEPL);
    fprintf('  perPercent = %.5g         graft load %.3f%% (target %.3f%%)\n', perPercent, out.load_wt, TGT_PCT);
    fprintf('  m_seed     = %.5g         BLZ945 load ratio %.3f (target %.3f)\n', m_seed, out.blz, TGT_BLZ);
    fprintf('  Q0les      = %.5g 1/day   lesion nu ratio %.3f (target %.2f)\n', Q0les, out.nuRatio, TGT_NU);
    fprintf('\n=== predictions (nothing below was fitted) ===\n');
    fprintf('  Irf8-/- graft cell count ratio : %.3f (measured %.3f)\n', out.nuRatio_ko, D.irf8Nu.rel);
    fprintf('  Irf8-/- graft load ratio       : %.3f (measured %.3f +- %.3f)\n', out.ratio, D.irf8Load.rel, D.irf8Load.relErr);
    fprintf('  host amyloid load at 4 weeks   : %.2f%% of area\n', out.hostPct);
    fprintf('  load-vs-carrier exponent gamma : %.3f (measured 0.687 +- 0.123)\n', out.gamma);
end
end

% =========================================================================
function p = base(variant, NC, k_a, pp, ms, cs, fn)
p = p210_params(variant);
p.N = NC;  p.k_a = k_a;  p.perPercent = pp;  p.m_seed = ms;
if ~isempty(cs), p.chi_s = cs; end
if ~isempty(fn), p.fnu_graft = fn; end
p = p210_derive(p);
end

function dep = depletion(NC, TS, k_a, pp, ms, cs, fn)
p   = base('baseline', NC, k_a, pp, ms, cs, fn);
sol = p210_solve(p, TS);  G = sol.G;
dep = 1 - sum(G.Vc(G.inHost) .* sol.m(end,G.inHost)') / G.Vhost;
end

function [L, H, nuG] = fullOf(variant, A)
p = base(variant, A.NC, A.k_a, A.pp, A.ms, A.cs, A.fn);
d = p210_post(p210_solve(p, A.TS));
L = d.loadPct(end);  H = d.hostPct(end);  nuG = d.nuGraft(end);
end

function L = loadOf(variant, NC, TS, k_a, pp, ms, cs, fn)
p = base(variant, NC, k_a, pp, ms, cs, fn);
d = p210_post(p210_solve(p, TS));
L = d.loadPct(end);
end

function nu = nuAt(day, NC, TS, k_a, pp, ms, cs, fn)
p = base('baseline', NC, k_a, pp, ms, cs, fn);
d = p210_post(p210_solve(p, TS));
nu = interp1(TS, d.nuGraft, day);
end

function R = blzRatio(NC, TS, k_a, pp, ms, cs, fn)
w = loadOf('baseline', NC, TS, k_a, pp, ms, cs, fn);
d = loadOf('depleted', NC, TS, k_a, pp, ms, cs, fn);
R = d / max(w, eps);
end

function R = agedRatio(NC, TS, k_a, pp, ms, cs, fn)
pw = base('baseline', NC, k_a, pp, ms, cs, fn);
pa = base('aged',     NC, k_a, pp, ms, cs, fn);
w = p210_post(p210_solve(pw, TS));
a = p210_post(p210_solve(pa, TS));
R = a.loadPct(end) / max(w.loadPct(end), eps);
end

function R = lesionRatio(Q0les, perPercent, chi_s, m_seed)
L = p210_lesion('baseline', Q0les, perPercent, chi_s, m_seed);
R = L.nuWin(end) / L.nuWin(1);
end

% -------------------------------------------------------------------------
function x = bisect(f, lo, hi, tol, itmax, verbose, name)
flo = f(lo);  fhi = f(hi);
% An endpoint that already meets the tolerance is the answer, whether or not
% the bracket changes sign across it.
if abs(flo) < tol
    if verbose, fprintf('  %s: lower end %g already within tolerance (%+.4g)\n', name, lo, flo); end
    x = lo;  return
end
if abs(fhi) < tol
    if verbose, fprintf('  %s: upper end %g already within tolerance (%+.4g)\n', name, hi, fhi); end
    x = hi;  return
end
if flo * fhi > 0
    if verbose
        fprintf('  %s: no sign change on [%g %g] (residuals %+.4g, %+.4g); keeping the closer end\n', ...
                name, lo, hi, flo, fhi);
    end
    if abs(flo) < abs(fhi), x = lo; else, x = hi; end
    return
end
x = sqrt(lo*hi);
for it = 1:itmax
    x  = sqrt(lo*hi);
    fx = f(x);
    if verbose, fprintf('  %-6s it %2d: %-11.5g residual %+.5g\n', name, it, x, fx); end
    if abs(fx) < tol, break; end
    if fx * flo > 0, lo = x; flo = fx; else, hi = x; end
end
end
