function R = p210_lesioncases()
%P210_LESIONCASES  The laser-lesion comparison with genotype-specific baselines.
%
%   R = P210_LESIONCASES
%
%   Fig. 3c of d'Errico et al. gives the amyloid area in the counting window
%   at 24 h, and Fig. 3d its increase over each genotype's own starting
%   level.  Together they imply different starting levels:
%       Irf8+/+   153.70 / 9.402 = 16.35 um^2 = 0.073% of the window
%       Irf8-/-    10.22 / 6.452 =  1.58 um^2 = 0.007% of the window
%   The runs of the paper start BOTH genotypes at 0.073% and compare fold
%   changes.  Dividing by the baseline does not remove that choice from a
%   nonlinear model, so this routine repeats the Irf8-/- run with its own
%   baseline and with the other things that differ in the knockout.
%
%   In every Irf8-/- case the window count at t = 0 is 0.825 of control and
%   the multiplier on chemotaxis is refitted so that the microglia in the
%   window rise 1.18-fold in 24 h, as measured (Fig. 3e).  Varied are
%     window baseline    0.073% (the paper), the Irf8-/- value 0.007%, and
%                        0.023%, their geometric mean, because the Irf8-/-
%                        value is a ratio of two means of which one is
%                        6.45 +/- 6.24
%     cargo per cell     the changes that bring it to 0.663 of control in
%                        the grafting model (p210_irf8): uptake divided by
%                        1.76, or k_d multiplied by 10, or k_r by 1.94
%     random motility    D_n multiplied by the same factor as chemotaxis
%     surrounding load   3.5% (the paper) or lower in proportion to the
%                        window baseline; ED Fig. 6b reports similar plaque
%                        counts in the two genotypes, so this last case is
%                        a what-if that the data do not support
%
%   The Irf8+/+ run is also repeated with its baseline halved and doubled.
%   Three numbers are compared with the measurement: the fold change of the
%   window amyloid at 24 h, the ratio of the fold changes of the two
%   genotypes, and the ratio of their ABSOLUTE amounts at 24 h, which does
%   not depend on any baseline estimate (measured 15.0).

D   = p210_data();
b0  = D.lesAb.t0 / D.lesAb.winUm2 * 100;          % % of area at t = 0: +/+, -/-
bm  = sqrt(b0(1) * b0(2));
p0  = p210_params('baseline');

fprintf('\n=== laser lesion: genotype-specific baselines ===\n');
fprintf('  window amyloid at t = 0 implied by Fig. 3c,d: Irf8+/+ %.4f%%, Irf8-/- %.4f%% of area\n', b0);
fprintf('  measured: fold change at 6 h and 24 h  Irf8+/+ %.2f +/- %.2f and %.2f +/- %.2f\n', ...
        D.lesInc.wt(2), D.lesInc.wtErr(2), D.lesInc.wt(3), D.lesInc.wtErr(3));
fprintf('                                         Irf8-/- %.2f +/- %.2f and %.2f +/- %.2f\n', ...
        D.lesInc.ko(2), D.lesInc.koErr(2), D.lesInc.ko(3), D.lesInc.koErr(3));
fprintf('            ratio of 24-h fold changes %.2f; ratio of 24-h areas %.1f; ratio of 24-h increments %.1f\n', ...
        D.lesAb.incRatio, D.lesAb.ratio, (D.lesAb.val(1) - D.lesAb.t0(1)) / (D.lesAb.val(2) - D.lesAb.t0(2)));

% ---- Irf8+/+ ---------------------------------------------------------------
fprintf('\n--- Irf8+/+ run: its own baseline, halved and doubled ---\n');
fprintf('%-28s %9s %9s %12s\n', 'window baseline, % area', 'fold 6 h', 'fold 24 h', 'amount 24 h');
Lw = [];
for f = [1 0.5 2]
    L  = p210_lesion('baseline', [], [], [], [], [], struct('winPct', f * b0(1)));
    i6 = find(L.tHr >= 6, 1, 'first');
    fprintf('%-28.4f %9.2f %9.2f %12.3f\n', f * b0(1), L.AbWinRel(i6), L.AbWinRel(end), L.AbWin(end));
    if f == 1, Lw = L; end
end
incW = Lw.AbWin(end) - Lw.AbWin(1);

% ---- Irf8-/- ---------------------------------------------------------------
C = { ...
  {'as in the paper (same baseline)',          struct('winPct', b0(1))}, ...
  {'own baseline',                             struct('winPct', b0(2))}, ...
  {'intermediate baseline',                    struct('winPct', bm)}, ...
  {'own baseline, lower uptake',               struct('winPct', b0(2), 'kuFac', 1.76)}, ...
  {'own baseline, faster degradation',         struct('winPct', b0(2), 'kdFac', 10)}, ...
  {'own baseline, faster release',             struct('winPct', b0(2), 'krFac', 1.94)}, ...
  {'own baseline, reduced motility',           struct('winPct', b0(2), 'fDn', NaN)}, ...
  {'own, reduced motility, lower uptake',      struct('winPct', b0(2), 'fDn', NaN, 'kuFac', 1.76)}, ...
  {'own, reduced motility, faster degr.',      struct('winPct', b0(2), 'fDn', NaN, 'kdFac', 10)}, ...
  {'own baseline, surrounding load 0.34%',     struct('winPct', b0(2), 'surPct', 3.5 * b0(2) / b0(1))}, ...
  {'own, surr. 0.34%, reduced motility',       struct('winPct', b0(2), 'surPct', 3.5 * b0(2) / b0(1), 'fDn', NaN)}};
fprintf('\n--- Irf8-/- run ---\n');
fprintf('%-38s %7s %7s | %8s %8s | %9s %9s %9s | %8s %8s\n', '', 'chemo-', 'cargo', 'fold', 'fold', 'ratio of', 'ratio of', 'ratio of', 'nu', 'nu');
fprintf('%-38s %7s %7s | %8s %8s | %9s %9s %9s | %8s %8s\n', '', 'taxis x', '/cell', '6 h', '24 h', 'folds', 'amounts', 'increm.', '6h/0h', '24h/0h');
R = struct('name',{},'fchi',{},'beta',{},'fold6',{},'fold24',{},'rFold',{},'rAmount',{},'rIncr',{},'nu6',{},'nu24',{});
for q = 1:numel(C)
    opt = C{q}{2};  tied = isfield(opt,'fDn');
    g = @(f) rise(f, opt, tied) - D.lesNu.koRel24;
    f = bis(g, 1e-3, 1.0, 5e-4);
    if tied, opt.fDn = f; end
    L  = p210_lesion('irf8ko', [], [], [], [], f, opt);
    i6 = find(L.tHr >= 6, 1, 'first');
    r.name = C{q}{1};  r.fchi = f;  r.beta = L.beta0 / Lw.beta0;
    r.fold6 = L.AbWinRel(i6);  r.fold24 = L.AbWinRel(end);
    r.rFold = Lw.AbWinRel(end) / L.AbWinRel(end);
    r.rAmount = Lw.AbWin(end) / L.AbWin(end);
    r.rIncr = incW / (L.AbWin(end) - L.AbWin(1));
    r.nu6 = L.nuWin(i6) / L.nuWin(1);  r.nu24 = L.nuWin(end) / L.nuWin(1);
    R(end+1) = r; %#ok<AGROW>
    fprintf('%-38s %7.3f %7.3f | %8.2f %8.2f | %9.2f %9.1f %9.1f | %8.2f %8.2f\n', r.name, r.fchi, r.beta, ...
            r.fold6, r.fold24, r.rFold, r.rAmount, r.rIncr, r.nu6, r.nu24);
end
fprintf('%-38s %7s %7s | %8.2f %8.2f | %9.2f %9.1f %9.1f | %8.2f %8.2f\n', 'measured', '', '0.663', D.lesInc.ko(2), D.lesInc.ko(3), ...
        D.lesAb.incRatio, D.lesAb.ratio, (D.lesAb.val(1) - D.lesAb.t0(1)) / (D.lesAb.val(2) - D.lesAb.t0(2)), ...
        D.lesNu.ko(2) / D.lesNu.ko(1), D.lesNu.koRel24);
fprintf('  "chemotaxis x" is the refitted multiplier on chi_s (and on D_n in the reduced-motility\n');
fprintf('  cases).  Cargo per cell is that of the surrounding cortex at t = 0, relative to the\n');
fprintf('  Irf8+/+ run.  Ratios are Irf8+/+ over Irf8-/-: of the 24-h fold changes, of the\n');
fprintf('  absolute window amounts at 24 h, and of the amounts gained in 24 h.\n');
end

% =========================================================================
function v = rise(f, opt, tied)
if tied, opt.fDn = f; end
L = p210_lesion('irf8ko', [], [], [], [], f, opt);
v = L.nuWin(end) / L.nuWin(1);
end

function x = bis(f, lo, hi, tol)
flo = f(lo);  fhi = f(hi);
if abs(flo) < tol, x = lo; return; end
if abs(fhi) < tol, x = hi; return; end
if flo * fhi > 0
    if abs(flo) < abs(fhi), x = lo; else, x = hi; end
    return
end
for it = 1:50
    x  = sqrt(lo * hi);  fx = f(x);
    if abs(fx) < tol, break; end
    if fx * flo > 0, lo = x; flo = fx; else, hi = x; end
end
end
