function R = p210_irf8()
%P210_IRF8  Irf8-/- recipients with host and graft constrained separately.
%
%   R = P210_IRF8
%
%   In the paper the Irf8-/- GRAFT run changes one number: the microglial
%   density is 0.133 of control everywhere, host and graft alike.  The 0.133
%   was measured inside the graft (Fig. 2g).  For the host, d'Errico et al.
%   give no cortical density, but two measurements bear on it and neither
%   is near 0.133:
%       0.825   microglia in the cortical imaging window of the lesion
%               experiment at t = 0, Irf8-/- against Irf8+/+ (Fig. 3e)
%       0.355   microglia associated with each cortical plaque (ED Fig. 6f)
%   and Horiuchi et al. (2012) report largely normal colonisation, 1.0.
%   The same figure also shows less amyloid inside each cortical microglial
%   cell, 0.663 of control (ED Fig. 6h).
%
%   This routine keeps what was measured in the graft and varies what was
%   not.  In every refitted case the graft density at 4 weeks is brought to
%   0.133 of the control run.  Varied are
%
%   HOST DENSITY         0.133 (the paper), 0.355, 0.825, 1.0
%
%   HOW THE GRAFT DENSITY COMES TO BE LOW
%     capacity           the graft's carrying capacity for microglia is
%                        reduced (fnu_graft refitted), chemotaxis unchanged:
%                        cells enter as in controls and are then removed
%     capacity, weak     the same, with chemotaxis multiplied by fchiLes,
%       chemotaxis       the factor fitted to the Irf8-/- lesion run (0.19)
%     motility           the cells move less: random motility D_n and
%                        chemotaxis chi_s are reduced by a common refitted
%                        factor, the graft capacity being that of controls
%
%   WHAT HAPPENS TO THE CARGO OF CELLS THAT ARE REMOVED
%     stays              as in the paper: the logistic term removes cells
%                        above the local capacity and leaves their cargo
%     leaves             p.cargoWithCells = 1: the cargo goes with them.
%                        k_a, perPercent and m_seed are refitted to their
%                        three targets for this variant, and its ratios are
%                        taken against its own control.
%
%   CARGO PER CELL       unchanged, or brought to 0.663 of control in the
%                        host at 4 weeks in one of three ways: lower uptake
%                        (k_u and phi k_p), faster degradation (k_d), or
%                        faster release (k_r).  Shown for two cases.
%
%   Nothing is fitted to the Irf8-/- graft LOAD in any case.

D = p210_data();
T = [0 14 28 112];
hosts = [0.355 0.825 1.0];

fprintf('\n=== Irf8-/-: host and graft constrained separately ===\n');
B = cell(1,2);
for cw = 0:1
    p = p210_params('baseline');
    if cw
        [ka, pp, ms] = refit(struct('cargoWithCells', 1), D);
        p.cargoWithCells = 1;  p.k_a = ka;  p.perPercent = pp;  p.m_seed = ms;  p = p210_derive(p);
    end
    s = p210_solve(p, T);  d = p210_post(s);
    B{cw+1} = struct('p', p, 'd', d, 'beta', hostBeta(s, 3));
end
p0 = B{1}.p;  dC = B{1}.d;
tgtNu = D.irf8Nu.rel * dC.nuGraft(3);
fprintf('  control: graft microglia %.0f and %.0f /mm^2 at 2 and 4 weeks, graft load %.3f%%, host load %.2f%%\n', ...
        dC.nuGraft(2), dC.nuGraft(3), dC.loadPct(3), dC.hostPct(3));
fprintf('  control, cargo leaves with removed cells: load %.3f%% before refitting; refitted k_a %.4f, perPercent %.2f, m_seed %.3g;\n', ...
        loadNoRefit(p0, T), B{2}.p.k_a, B{2}.p.perPercent, B{2}.p.m_seed);
fprintf('    load %.3f%%, host load %.2f%%, carrier share %.1f%% at 4 weeks and %.1f%% at 16 weeks\n', ...
        B{2}.d.loadPct(3), B{2}.d.hostPct(3), 100*B{2}.d.carrShare(3), 100*B{2}.d.carrShare(4));
fprintf('  target graft density at 4 weeks: %.3f of control = %.0f /mm^2 (Fig. 2g)\n', D.irf8Nu.rel, tgtNu);
fprintf('  measured graft load ratio: %.3f +/- %.3f (Fig. 2h)\n', D.irf8Load.rel, D.irf8Load.relErr);

% ---- the fits of the microglial field (they do not depend on the amyloid)
fit = struct('h',{},'entry',{},'x',{},'dn',{},'fn',{},'ok',{});
for h = hosts
    for x = [1 p0.fchiLes]
        g = @(f) nu4(setIrf8(p0, h, x, 1, f, 'unchanged', 1), T) - tgtNu;
        [fn, ok] = bis(g, 1e-4, 3, 0.2);
        if x == 1, nm = 'capacity'; else, nm = 'capacity, weak chemotaxis'; end
        fit(end+1) = struct('h',h,'entry',nm,'x',x,'dn',1,'fn',fn,'ok',ok); %#ok<AGROW>
    end
    g = @(a) nu4(setIrf8(p0, h, a, a, p0.fnu_graft, 'unchanged', 1), T) - tgtNu;
    [a, ok] = bis(g, 1e-5, 1, 0.2);
    fit(end+1) = struct('h',h,'entry','motility','x',a,'dn',a,'fn',p0.fnu_graft,'ok',ok); %#ok<AGROW>
end

R = struct('host',{},'entry',{},'cargoFate',{},'cargo',{},'par',{},'nu2',{},'nu4',{},'ratio',{},'import',{},'hostRel',{},'betaRel',{});
hdr = @() fprintf('\n%-6s %-26s %-7s %-19s | %9s %7s %7s | %7s %8s %7s %9s\n%-6s %-26s %-7s %-19s | %9s %7s %7s | %7s %8s %7s %9s\n', ...
      'host', 'graft density low by', 'cargo', 'cargo per cell', 'refitted', 'nu 2wk', 'nu 4wk', 'LOAD', 'carrier', 'host', 'cargo', ...
      'dens.', '', 'of lost', '', 'factor', '/mm^2', '/mm^2', 'RATIO', 'import', 'load', 'per cell');
fate = {'stays', 'leaves'};

fprintf('\n--- A. host density and the way the graft density comes to be low ---');
hdr();
for cw = 0:1
    b = B{cw+1};
    pk = b.p;  pk.nu0 = pk.nu0IRF8;  pk = p210_derive(pk);          % the reading of the paper
    R(end+1) = report(pk, b, T, pk.nu0IRF8, 'density 0.133 everywhere', fate{cw+1}, 'unchanged', NaN); %#ok<AGROW>
    for k = 1:numel(fit)
        f = fit(k);
        pv = setIrf8(b.p, f.h, f.x, f.dn, f.fn, 'unchanged', 1);
        if strcmp(f.entry, 'motility'), par = f.x; else, par = f.fn; end
        R(end+1) = report(pv, b, T, f.h, f.entry, fate{cw+1}, 'unchanged', par); %#ok<AGROW>
        if ~f.ok, fprintf('   (graft density target not attainable)\n'); end
    end
end
fprintf('  "refitted factor" is the graft capacity fnu_graft (control %.3f) in the capacity cases\n', p0.fnu_graft);
fprintf('  and the common multiplier on D_n and chi_s in the motility cases.  LOAD RATIO is the\n');
fprintf('  graft load at 4 weeks over control (measured %.3f +/- %.3f); carrier import is the\n', D.irf8Load.rel, D.irf8Load.relErr);
fprintf('  cumulative carrier-borne import into the graft by 4 weeks over control; host load and\n');
fprintf('  cargo per cell in the host are also relative to control.\n');

fprintf('\n--- B. cargo per cell brought to %.3f of control (ED Fig. 6h), three ways ---', p0.cargoIRF8);
hdr();
cargos = {'lower uptake', 'faster degradation', 'faster release'};
k825 = find([fit.h] == 0.825 & strcmp({fit.entry}, 'capacity, weak chemotaxis'), 1);
cases = {{1, p0.nu0IRF8, 1, 1, p0.fnu_graft, 'density 0.133 everywhere'}, ...
         {2, fit(k825).h, fit(k825).x, 1, fit(k825).fn, fit(k825).entry}};
for q = 1:numel(cases)
    cs = cases{q};  b = B{cs{1}};
    for j = 1:numel(cargos)
        gb = @(a) hostBeta(p210_solve(setIrf8(b.p, cs{2}, cs{3}, cs{4}, cs{5}, cargos{j}, a), T), 3) / b.beta - p0.cargoIRF8;
        [a, ok] = bis(gb, 1, 60, 2e-3);
        pv = setIrf8(b.p, cs{2}, cs{3}, cs{4}, cs{5}, cargos{j}, a);
        R(end+1) = report(pv, b, T, cs{2}, cs{6}, fate{cs{1}}, cargos{j}, a); %#ok<AGROW>
        if ~ok, fprintf('   (cargo target not attainable)\n'); end
    end
end
fprintf('  here "refitted factor" is the factor by which the uptake is divided, or k_d or k_r\n');
fprintf('  multiplied.\n');

r = [R.ratio];
fprintf('\n  Load ratio over all cases: %.3f to %.3f; measured %.3f +/- %.3f.\n', min(r), max(r), D.irf8Load.rel, D.irf8Load.relErr);
im = [R.import];
cc = corrcoef(log(im), log(r));
fprintf('  Across the cases the load ratio follows the carrier-borne import: correlation of the\n');
fprintf('  logarithms %.2f, while the graft density at 4 weeks is the same in all of them.\n', cc(1,2));
end

% =========================================================================
function p = setIrf8(pb, h, x, dn, fn, cargo, a)
p = pb;
p.nu0 = h;  p.chi_s = x * pb.chi_s;  p.D_n = dn * pb.D_n;  p.fnu_graft = fn;
switch cargo
    case 'unchanged'
    case 'lower uptake',       p.k_u = pb.k_u / a;  p.phi = pb.phi / a;
    case 'faster degradation', p.k_d = pb.k_d * a;
    case 'faster release',     p.k_r = pb.k_r * a;
end
p = p210_derive(p);
end

function v = nu4(p, T)
d = p210_post(p210_solve(p, T));
v = d.nuGraft(3);
end

function L = loadNoRefit(p0, T)
p = p0;  p.cargoWithCells = 1;
d = p210_post(p210_solve(p, T));
L = d.loadPct(3);
end

function b = hostBeta(sol, k)
G = sol.G;  hi = G.inHost;
b = sum(G.Vc(hi) .* sol.b(k,hi)') / sum(G.Vc(hi) .* sol.nu(k,hi)');
end

function r = report(p, b, T, h, entry, fate, cargo, par)
s = p210_solve(p, T);  d = p210_post(s);
r.host = h;  r.entry = entry;  r.cargoFate = fate;  r.cargo = cargo;  r.par = par;
r.nu2 = d.nuGraft(2);  r.nu4 = d.nuGraft(3);
r.ratio   = d.loadPct(3) / b.d.loadPct(3);
r.import  = d.carrImport(3) / b.d.carrImport(3);
r.hostRel = d.hostPct(3) / b.d.hostPct(3);
r.betaRel = hostBeta(s, 3) / b.beta;
fprintf('%-6.3f %-26s %-7s %-19s | %9.4f %7.0f %7.0f | %7.3f %8.3f %7.3f %9.3f\n', ...
        h, entry, fate, cargo, par, r.nu2, r.nu4, r.ratio, r.import, r.hostRel, r.betaRel);
end

function [ka, pp, ms] = refit(ov, D)
%Refit k_a, perPercent and m_seed to their three targets, as p210_calibrate
%does, for a variant of the model given by the fields of ov.
NC = 400;  TS = [0 7 14 21 28];
L0 = load('p210_fitted.mat');
ka = L0.k_a;  pp = L0.perPercent;  ms = L0.m_seed;
for outer = 1:16
    fa = @(k) depl(mk('baseline',ov,k,pp,ms,NC), TS) - 0.80;
    ka = bis(fa, 1e-3, 5, 1e-4);
    ppOld = pp;
    for it = 1:40
        L  = loadOf(mk('baseline',ov,ka,pp,ms,NC), TS);
        np = pp * (L / D.irf8Load.val(1))^0.75;  rel = abs(np - pp) / np;  pp = np;
        if rel < 5e-5, break; end
    end
    fm = @(m) loadOf(mk('depleted',ov,ka,pp,m,NC), TS) / max(loadOf(mk('baseline',ov,ka,pp,m,NC), TS), eps) - D.deplete.rel;
    msNew = bis(fm, 1e-5, 3e-3, 1e-4);
    relm  = abs(msNew - ms) / msNew;  ms = msNew;
    if abs(pp - ppOld)/pp < 2e-4 && relm < 5e-4, break; end
end
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
    ok = false;
    if abs(flo) < abs(fhi), x = lo; else, x = hi; end
    return
end
for it = 1:50
    x  = sqrt(lo * hi);  fx = f(x);
    if abs(fx) < tol, break; end
    if fx * flo > 0, lo = x; flo = fx; else, hi = x; end
end
end
