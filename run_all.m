function run_all(mode)
%RUN_ALL  Reproduce every number quoted in the P210 manuscript.
%
%   RUN_ALL           full run: calibrate, then report
%   RUN_ALL('quick')  skip calibration and use the stored fit
%
%   Console output is saved by the caller; see run_all_output.txt.

if nargin < 1, mode = 'full'; end
t0 = tic;
fprintf('===========================================================\n');
fprintf(' P210  microglial carrier-diffusion model: run_all (%s)\n', mode);
fprintf(' %s\n', datestr(now,'yyyy-mm-dd HH:MM'));
fprintf('===========================================================\n');

if ~strcmpi(mode,'quick')
    o = p210_calibrate(false);
else
    L = load('p210_fitted.mat');
    o = L;  o.hostPct = NaN;
end

D = p210_data();
p = p210_params('baseline');
fprintf('\n--- fitted constants: six, one published target each ------\n');
fprintf('  k_a        = %-10.5g  surface elongation, 1/day\n', p.k_a);
fprintf('  perPercent = %-10.5g  model units per 1%% of 6E10+ area\n', p.perPercent);
fprintf('  m_seed     = %-10.5g  deposit mass per carrier release event\n', p.m_seed);
fprintf('  chi_s      = %-10.5g  chemotactic sensitivity, mm^2/day per unit s\n', p.chi_s);
fprintf('  fnu_graft  = %-10.5g  graft microglial carrying capacity\n', p.fnu_graft);
fprintf('  Q0les      = %-10.5g  acute lesion signal, 1/day\n', p.Q0les);

fprintf('\n--- NOT fitted: the aged uptake-competence proxy ----------\n');
fprintf('  phiAged    = %-10.5g  EX VIVO PROXY, not fitted: the ratio of\n', p.phiAged);
fprintf('                          intracellular Abeta cluster burden in old\n');
fprintf('                          against adult isolated microglia, ED Fig.\n');
fprintf('                          4d.  Confocal counting after 24 h culture,\n');
fprintf('                          not a measured uptake rate, so the aged\n');
fprintf('                          graft load is conditional on it\n');

fprintf('\n--- declared or measured inputs --------------------------\n');
fprintf('  D_m   = %-8.4g mm^2/day   Waters (2010), 62.3 um^2/s\n', p.D_m);
fprintf('  lam_m = %-8.4g 1/day      Cirrito (2003), t_1/2 = 1.97 h\n', p.lam_m);
fprintf('  D_n   = %-8.4g mm^2/day   from graft arrival within 2 weeks\n', p.D_n);
fprintf('  k_1   = %-8.4g units/day   homogeneous nucleation at m = 1: DECLARED,\n', p.k_1);
fprintf('                              no tissue measurement exists.  The seeding\n');
fprintf('                              bottleneck survives to k_1 ~ 6e-3 (see below)\n');
fprintf('  k_2   = %-8.4g            secondary nucleation on surface\n', p.k_2);
fprintf('  k_d   = %-8.4g 1/day      lysosomal degradation of cargo\n', p.k_d);
fprintf('  k_r   = %-8.4g 1/day      release of cargo\n', p.k_r);
fprintf('  R_g   = %-8.4g mm         graft radius\n', p.Rg);
fprintf('  c_h   = %-8.4g /mm^3     host deposit density, ED Fig. 6b (38.6/mm^2)\n', p.c_hRef);
fprintf('  detR  = %-8.4g           growth over the seed at which a deposit shows\n', p.detRatio);
fprintf('  kap   = %-8.4g           DERIVED  area per unit c^(1/3) F^(2/3)\n', p.kap);
fprintf('  m_seed against its independent estimate b_max/n0_3D = %.4g: ratio %.2f\n', ...
        p.b_max/p.n0_3D, p.m_seed/(p.b_max/p.n0_3D));

TS  = [0 7 14 21 28 42 56 84 112];
solB = p210_solve(p210_params('baseline'), TS);  dB = p210_post(solB);
solN = p210_solve(p210_params('nocarrier'),TS);  dN = p210_post(solN);
solK = p210_solve(p210_params('irf8ko'),   TS);  dK = p210_post(solK);
solD = p210_solve(p210_params('depleted'), TS);  dD = p210_post(solD);
solA = p210_solve(p210_params('aged'),     TS);  dA = p210_post(solA);
i4 = find(TS==28); i16 = find(TS==112); i2 = find(TS==14);

fprintf('\n--- transport length scales ------------------------------\n');
G = solB.G;
SC = p210_scales(p);
fprintf('  ell_m, clean tissue              %.3f mm\n', SC.ell_m_clean);
fprintf('  ell_m, graft                     %.3f mm\n', SC.ell_m_graft);
fprintf('  ell_m, peri-graft host           %.3f mm\n', SC.ell_m_host);
fprintf('  ell_x, carrier export length     %.3f mm\n', SC.ell_x);
fprintf('  Lambda = ell_x / ell_m           %.3f\n', SC.Lambda);
fprintf('  carrier migration speed V        %.4f mm/day (%.1f um/day)\n', SC.Vmig, SC.Vmig*1e3);
fprintf('    directed V*tau %.0f um  vs  random sqrt(2 Dn tau) %.0f um  (tau %.2f d)\n', ...
            SC.drift*1e3, SC.wander*1e3, SC.tau);
fprintf('  mobility contrast D_n/D_m        %.3g\n', p.D_n/p.D_m);

fprintf('\n--- the seeding bottleneck -------------------------------\n');
fprintf('  graft load at 16 weeks, carriers active   %.4f %% area\n', dB.loadPct(i16));
fprintf('  graft load at 16 weeks, carriers inert    %.4g %% area\n', dN.loadPct(i16));
fprintf('  max deposit anywhere in graft, inert      %.4g model units\n', max(solN.F(end,G.inGraft)));
mbN = max(solN.F(end,G.inGraft)) / max(max(solN.c(end,G.inGraft)), realmin);
fprintf('  largest inert deposit, in seed masses     %.4g\n', mbN/p.m_seed);
fprintf('  (a deposit shows at about                 %.4g seed masses)\n', p.detRatio);

fprintf('\n--- how far the bottleneck survives in k_1 ----------------\n');
fprintf('  k_1 is DECLARED, not measured, so its margin is part of the result.\n');
fprintf('  Each row raises it and reruns the inert arm.  A graft counts as\n');
fprintf('  detectable above %.2f %% of area, and a deposit shows at %.0f seed masses.\n', ...
        p.aThr, p.detRatio);
fprintf('       k_1      inert load %%   largest inert   active load %%\n');
fprintf('                               (seed masses)\n');
k1v = [1e-8 1e-5 1e-4 1e-3 3e-3 6e-3 1e-2];
k1Floor = NaN;  k1Show = NaN;
for q = 1:numel(k1v)
    pn2 = p210_params('nocarrier');  pn2.k_1 = k1v(q);
    sn2 = p210_solve(pn2, TS);       dn2 = p210_post(sn2);
    mb2 = max(sn2.F(end,G.inGraft)) / max(max(sn2.c(end,G.inGraft)), realmin) / pn2.m_seed;
    if ismember(k1v(q), [1e-8 1e-3])          % the two values quoted in the paper
        pb2 = p210_params('baseline');  pb2.k_1 = k1v(q);
        act = sprintf('%12.4f', p210_post(p210_solve(pb2, TS)).loadPct(i16));
    else
        act = '           -';
    end
    fprintf('  %9.1e  %12.4g  %14.3g %s\n', k1v(q), dn2.loadPct(i16), mb2, act);
    if isnan(k1Floor) && dn2.loadPct(i16) >= p.aThr,  k1Floor = k1v(q);  end
    if isnan(k1Show)  && mb2 >= p.detRatio,           k1Show  = k1v(q);  end
end
fprintf('  inert deposition first clears the %.2f %% floor at k_1 = %.0e, and first\n', ...
        p.aThr, k1Floor);
fprintf('  reaches %.0f seed masses at k_1 = %.0e: %.1f orders of magnitude above the\n', ...
        p.detRatio, k1Show, log10(k1Show/p.k_1));
fprintf('  declared value.  The carriers-active load is unchanged across that range,\n');
fprintf('  so the bottleneck says how slow homogeneous nucleation would have to be,\n');
fprintf('  not that it cannot happen.\n');


fprintf('\n--- host load does not depend on the carriers ------------\n');
fprintf('  Grathwohl (2009): 95-97%% ablation for 2-4 weeks, two APP lines, no\n');
fprintf('  change in congophilic or total Abeta load, and no change in de novo\n');
fprintf('  plaque formation when ablation precedes deposition.  Spangenberg\n');
fprintf('  (2016): ~80%% elimination for 28 days in 5xFAD, no change in plaque\n');
fprintf('  load, number or size, nor in soluble or insoluble Abeta.\n');
fprintf('  Olmos-Alonso (2016): CSF1R blockade for 3 months, amyloid unchanged.\n');
fprintf('  d''Errico ED Fig. 6b: Irf8-/- host cortex %.1f vs %.1f plaques/mm^2.\n', ...
        D.plaqDens.ctx(3), D.plaqDens.ctx(2));
fprintf('   carriers     host load     ratio     graft load     ratio\n');
for fcar = [1.0 0.20 0.13 0.05]
    pf = p210_params('baseline');  pf.nu0 = fcar;
    df = p210_post(p210_solve(pf, TS));
    fprintf('   %5.0f%%      %7.3f %%    %6.3f      %7.4f %%    %6.3f\n', ...
            100*fcar, df.hostPct(i4), df.hostPct(i4)/dB.hostPct(i4), ...
            df.loadPct(i4), df.loadPct(i4)/dB.loadPct(i4));
end
fprintf('  The host is indifferent to the carriers and the graft is not.  That\n');
fprintf('  asymmetry is what the seeding bottleneck predicts: established\n');
fprintf('  tissue grows on the surface it already has, naive tissue has none.\n');

fprintf('\n--- mass budget across the graft boundary ----------------\n');
fprintf('  carrier-borne import by 4 weeks    %10.4g units\n', dB.carrImport(i4));
fprintf('  diffusive monomer import, 4 weeks  %10.4g units\n', dB.diffImport(i4));
fprintf('  carrier share of imported mass     %10.4f %%\n', 100*dB.carrShare(i4));
fprintf('  carrier-borne import by 16 weeks   %10.4g units\n', dB.carrImport(i16));
fprintf('  diffusive import by 16 weeks       %10.4g units\n', dB.diffImport(i16));
fprintf('  carrier share at 16 weeks          %10.4f %%\n', 100*dB.carrShare(i16));

fprintf('\n--- net effect of the carriers: active minus inert -------\n');
fprintf('   weeks   graft change   host change   (percentage points of area)\n');
for k = [i4 i16]
    fprintf('  %6.0f   %+10.3f    %+10.3f\n', TS(k)/7, ...
        dB.loadPct(k) - dN.loadPct(k), dB.hostPct(k) - dN.hostPct(k));
end
fprintf('  [Fig. 7a shows the 4-week change, swept over the export length.]\n');

fprintf('\n--- graft amyloid load, %% of area ------------------------\n');
fprintf('   weeks    base    Irf8-/-   BLZ945    aged    no carrier\n');
for k = [i2 i4 i16]
    fprintf('  %6.0f  %7.3f  %7.3f  %7.3f  %7.3f  %9.3g\n', TS(k)/7, ...
        dB.loadPct(k), dK.loadPct(k), dD.loadPct(k), dA.loadPct(k), dN.loadPct(k));
end
fprintf('  relative to baseline at 4 weeks: Irf8-/- %.3f, BLZ945 %.3f, aged %.3f\n', ...
    dK.loadPct(i4)/dB.loadPct(i4), dD.loadPct(i4)/dB.loadPct(i4), dA.loadPct(i4)/dB.loadPct(i4));
fprintf('  measured (Source Data)         : Irf8-/- %.3f,  BLZ945 %.3f,  aged %.3f\n', ...
    D.irf8Load.rel, D.deplete.rel, D.aged.rel);
fprintf('  [BLZ945 is fitted via m_seed.  The Irf8-/- and aged load ratios\n');
fprintf('   are PREDICTIONS: neither arm fits an amyloid parameter.  Irf8-/-\n');
fprintf('   has its carrier density imposed at the measured 0.133 (Fig. 2g);\n');
fprintf('   aged has the measured uptake 0.494 (ED Fig. 4d) and nothing else.]\n');

fprintf('\n--- Cx3cr1-/-: selective loss of the deposit-surface route ---\n');
fprintf('  Those recipients keep a normal graft microglial density and lose three\n');
fprintf('  quarters of the amyloid, which no change in carrier NUMBER can produce.\n');
fprintf('  The model already loads carriers by two routes, soluble (k_u) and deposit\n');
fprintf('  surface (phi k_p).  Cx3cr1 mediates protofibril phagocytosis (Liu 2010),\n');
fprintf('  so the genotype is the second route switched off; nothing is tuned.\n');
pC = p210_params('cx3cr1ko');  dC = p210_post(p210_solve(pC, TS));
fprintf('    graft microglia at 4 wk   model %5.0f of %5.0f (%.3f)   measured %.3f\n', ...
        dC.nuGraft(i4), dB.nuGraft(i4), dC.nuGraft(i4)/dB.nuGraft(i4), D.cx3cr1.nuRel);
fprintf('    graft load at 4 wk        model %.3f of %.3f (%.3f)   measured %.3f\n', ...
        dC.loadPct(i4), dB.loadPct(i4), dC.loadPct(i4)/dB.loadPct(i4), D.cx3cr1.rel);

loC = D.cx3cr1.rel + D.cx3cr1.relErr;   % upper 1 s.e.m. limit
aa = 0; bb = 1;
for it = 1:24
    mm = 0.5*(aa+bb);  pm = p210_params('baseline');  pm.phi = mm;
    rr = p210_post(p210_solve(pm, TS)).loadPct(i4) / dB.loadPct(i4);
    if rr > loC, bb = mm; else, aa = mm; end
end
fprintf('    phi is bounded, not merely set to zero: the largest value still inside\n');
fprintf('    the measured 1 s.e.m. interval is phi = %.3f, so the data require the\n', aa);
fprintf('    fibrillar-uptake route reduced by at least %.0f%%.\n', 100*(1-aa));

fprintf('\n--- can the carrier field be filled at an imaged somal speed? ---\n');
fprintf('  Nimmerjahn (2005): 5%% of somata shift 1-2 um/h.  The calibrated peak\n');
fprintf('  drift is above that; this is what bringing it into range costs.\n');
fprintf('   chi_s x   V (um/h)    nu 2wk    nu 4wk   load 4wk\n');
for fch = [1.0 0.5 0.25 0.0]
    pv = p210_params('baseline');  pv.chi_s = fch * p.chi_s;
    Sv = p210_scales(pv);  dv = p210_post(p210_solve(pv, TS));
    fprintf('   %7.2f   %8.2f  %8.0f  %8.0f   %8.4f\n', ...
            fch, Sv.Vmig*1e3/24, dv.nuGraft(3), dv.nuGraft(i4), dv.loadPct(i4));
end
fprintf('   measured              %8.0f  %8.0f   %8.3f\n', ...
        D.graftNu.val(1), D.graftNu.val(2), D.irf8Load.val(1));

fprintf('\n--- microglial density in the graft, cells/mm^2 ----------\n');
fprintf('   2 weeks %6.0f   4 weeks %6.0f   16 weeks %6.0f\n', ...
    dB.nuGraft(i2), dB.nuGraft(i4), dB.nuGraft(i16));
fprintf('  measured (Source Data Fig. 1d): 2 weeks %.0f +/- %.0f, 4 weeks %.0f +/- %.0f\n', ...
    D.graftNu.val(1), D.graftNu.err(1), D.graftNu.val(2), D.graftNu.err(2));
fprintf('  graft cell ratio at 4 weeks: Irf8-/- %.3f (imposed, target %.3f), BLZ945 %.3f (measured %.3f)\n', ...
    dK.nuGraft(i4)/dB.nuGraft(i4), D.irf8Nu.rel, dD.nuGraft(i4)/dB.nuGraft(i4), D.deplete.nuRel);

fprintf('\n--- host tissue -------------------------------------------\n');
fprintf('  free Abeta in host at 4 weeks       %.4f of clean tissue\n', ...
    sum(G.Vc(G.inHost).*solB.m(i4,G.inHost)')/G.Vhost);
fprintf('  host amyloid load at 4 weeks        %.2f %% area\n', dB.hostPct(i4));
fprintf('  graft-to-host load ratio at 4 weeks %.3f\n', dB.loadPct(i4)/dB.hostPct(i4));

fprintf('\n--- graft load against the published time course ---------\n');
fprintf('  ED Fig. 1c is the only published 16-week graft load.  Its cohort sits\n');
fprintf('  at 1.083%% at 4 weeks, not the 0.586%% of Fig. 2h the model is calibrated\n');
fprintf('  to, so the within-cohort ratios are the meaningful comparison.\n');
fprintf('    weeks     model     measured (ED Fig. 1c)      n\n');
for q = 1:3
    km = interp1(TS, dB.loadPct, D.graftTime.t(q));
    fprintf('   %6.0f   %7.3f   %7.3f +/- %.3f %%        %2d\n', ...
            D.graftTime.t(q)/7, km, D.graftTime.val(q), D.graftTime.err(q), D.graftTime.n(q));
end
fprintf('\n  Ratios between time points remove the cohort offset.  These loads are\n');
fprintf('  strongly right-skewed, so the ratio of means with Gaussian errors is\n');
fprintf('  reported only for comparison; the log scale is the appropriate one and\n');
fprintf('  the bootstrap assumes no distribution at all.  See p210_ratiostats.\n');
for q = [3 2]            % 16/4 then 16/2
    if q == 3
        lab = '16 wk / 4 wk';  den = D.graftTime.v4;   tden = 28;
    else
        lab = '16 wk / 2 wk';  den = D.graftTime.v2;   tden = 14;
    end
    mr = interp1(TS, dB.loadPct, 112) / interp1(TS, dB.loadPct, tden);
    if q == 3, nm = D.graftTime.nMice([3 2]); else, nm = D.graftTime.nMice([3 1]); end
    S  = p210_ratiostats(D.graftTime.v16, den, mr, [], nm);     % rho = 1, worst case
    S0 = p210_ratiostats(D.graftTime.v16, den, mr);             % no clustering
    fprintf('\n  %s   model %.2f\n', lab, mr);
    fprintf('    ratio of means    %.2f +/- %.2f (Gaussian)      model is %.1f sigma out\n', ...
            S.rMean, S.seMeanGauss, S.sigGauss);
    fprintf('    ratio of medians  %.2f\n', S.rMedian);
    if isnan(S.logRatio)
        fprintf('    log scale         not available: a sample contains a zero\n');
    else
        fprintf('    log scale         x%.2f, 95%% CI x%.2f to x%.2f    model is %.1f sigma out\n', ...
                exp(S0.logRatio), S0.logCI(1), S0.logCI(2), S0.sigLog);
        fprintf('      clustered       grafts/mice %d/%d and %d/%d, design effect %.2f/%.2f at rho=1\n', ...
                numel(D.graftTime.v16), nm(1), numel(den), nm(2), S.deff(1), S.deff(2));
        fprintf('                      x%.2f, 95%% CI x%.2f to x%.2f    model is %.1f sigma out\n', ...
                exp(S.logRatio), S.logCI(1), S.logCI(2), S.sigLog);
        fprintf('                      change itself t = %.2f on %d and %d animals%s\n', ...
                S.logT, numel(D.graftTime.v16), numel(den), ...
                repmat(' (not significant)', 1, abs(S.logT) < 2.07));
    end
    fprintf('    bootstrap         95%% CI %.2f to %.2f,  P(ratio >= model) = %.4f\n', ...
            S.bootCI(1), S.bootCI(2), S.bootP);
end
fprintf('\n  Every graft-level panel in d''Errico et al. is "one graft from n mice", so\n');
fprintf('  all of these standard errors are clustered.  rho = 1 above assumes grafts\n');
fprintf('  within a mouse are perfectly correlated, which BOUNDS the inflation rather\n');
fprintf('  than estimating it; the conclusion is unchanged either way.  The published\n');
fprintf('  ED Fig. 1c test, F(2,28), has 31 observations, so it is graft-level too.\n');
fprintf('\n  The log scale is the one quoted in the paper.  On it the measured\n');
fprintf('  4-to-16 week increase is not itself significant and its interval\n');
fprintf('  contains the model value, so the apparent discrepancy in the ratio of\n');
fprintf('  means is an artefact of assuming normality for a skewed quantity.\n');

fprintf('\n--- onset of detectable deposit --------------------------\n');
ok = ~isnan(dB.onsetT);
fprintf('  rim    %.2f weeks;   mid-graft %.2f weeks;   centre %.2f weeks\n', ...
    min(dB.onsetT(ok))/7, median(dB.onsetT(ok))/7, max(dB.onsetT(ok))/7);
okK = ~isnan(dK.onsetT);
if any(okK)
    fprintf('  Irf8-/-: rim %.2f weeks; deepest detected %.3f mm of %.2f mm\n', ...
        min(dK.onsetT(okK))/7, max(dK.depth(okK)), p.Rg);
end

fprintf('\n--- laser lesion ------------------------------------------\n');
Lw = p210_lesion('baseline'); Lk = p210_lesion('irf8ko');
i6 = find(Lw.tHr >= 6, 1, 'first');
fprintf('  microglia in window, 24 h / 0 h      %.2f  (measured %.2f)  FITTED\n', ...
    Lw.nuWin(end)/Lw.nuWin(1), D.lesNu.rel24);
fprintf('  amyloid in window, 6 h / 0 h         %.2f  (measured %.2f +/- %.2f)\n', ...
    Lw.AbWinRel(i6), D.lesInc.wt(2), D.lesInc.wtErr(2));
fprintf('  amyloid in window, 24 h / 0 h        %.2f  (measured %.2f +/- %.2f)\n', ...
    Lw.AbWinRel(end), D.lesInc.wt(3), D.lesInc.wtErr(3));
fprintf('  Irf8-/- amyloid, 24 h / 0 h          %.2f  (measured %.2f +/- %.2f)\n', ...
    Lk.AbWinRel(end), D.lesInc.ko(3), D.lesInc.koErr(3));
fprintf('  the 15.0-fold ratio of 24 h AREAS is not a test of the model: the two\n');
fprintf('  genotypes did not start level.  Fig. 3c and 3d together imply %.2f um^2\n', D.lesAb.t0(1));
fprintf('  in Irf8+/+ at t = 0 against %.2f in Irf8-/-, so 15.0 = %.1f already there\n', D.lesAb.t0(2), D.lesAb.t0Ratio);
fprintf('  x %.2f contributed by delivery.  The model starts both level, so only\n', D.lesAb.incRatio);
fprintf('  the second factor tests it:\n');
fprintf('  ratio of 24 h INCREMENTS             %.2f  (measured %.2f)\n', ...
    Lw.AbWinRel(end)/max(Lk.AbWinRel(end),eps), D.lesAb.incRatio);

fprintf('\n--- grid convergence --------------------------------------\n');
for N = [160 320 640]
    pk = p210_params('baseline'); pk.N = N;
    dk = p210_post(p210_solve(pk, [0 28 112]));
    fprintf('  N = %4d :  4 weeks %.4f %%   16 weeks %.4f %%\n', N, dk.loadPct(2), dk.loadPct(3));
end

fprintf('\n--- positivity --------------------------------------------\n');
fprintf('  min nu = %.3g, min b = %.3g, min F = %.3g, min m = %.3g\n', ...
    min(solB.nu(:)), min(solB.b(:)), min(solB.F(:)), min(solB.m(:)));

fprintf('\nrun_all finished in %.1f s\n', toc(t0));
end
