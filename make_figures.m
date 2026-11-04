function R = make_figures(outdir)
%MAKE_FIGURES  Produce Figs 1-8 and Figs S1-S3 of the P210 manuscript.
%
%   R = MAKE_FIGURES            writes to ../figures
%   R = MAKE_FIGURES(outdir)    writes to outdir
%
%   Returns a struct of the numbers quoted in the text, so that run_all can
%   print them without recomputing.

if nargin < 1 || isempty(outdir), outdir = fullfile('..','figures'); end
if ~exist(outdir,'dir'), mkdir(outdir); end
set(0,'DefaultAxesFontSize',9,'DefaultAxesFontName','Helvetica', ...
      'DefaultLineLineWidth',1.3);
D  = p210_data();
TS = [0 3 7 14 21 28 42 56 84 112];
R  = struct();

% =====================================================================
% baseline and perturbation runs
% =====================================================================
fprintf('running baseline and perturbations ...\n'); tic
solB = p210_solve(p210_params('baseline'),  TS);   dB = p210_post(solB);
solK = p210_solve(p210_params('irf8ko'),    TS);   dK = p210_post(solK);
solD = p210_solve(p210_params('depleted'),  TS);   dD = p210_post(solD);
solA = p210_solve(p210_params('aged'),      TS);   dA = p210_post(solA);
solN = p210_solve(p210_params('nocarrier'), TS);   dN = p210_post(solN);
fprintf('  done in %.1f s\n', toc);

p = solB.p;  G = solB.G;
i4  = find(TS==28);  i16 = find(TS==112);  i2 = find(TS==14);

R.load4   = dB.loadPct(i4);    R.load16 = dB.loadPct(i16);
R.loadKO4 = dK.loadPct(i4);    R.loadD4 = dD.loadPct(i4);
R.loadA4  = dA.loadPct(i4);    R.loadN16 = dN.loadPct(i16);
R.nu2     = dB.nuGraft(i2);    R.nu4    = dB.nuGraft(i4);
R.Lambda  = dB.Lambda;  R.ell_x = dB.ell_x;  R.Vmig = dB.Vmig;
SC = p210_scales(p);          % resolved on a fine time grid; see p210_scales
R.ell_m_graft = SC.ell_m_graft;
R.ell_m_host  = SC.ell_m_host;
R.ell_x = SC.ell_x;  R.Lambda = SC.Lambda;  R.Vmig = SC.Vmig;
R.carrFrac4   = dB.carrFrac(i4);  R.carrFrac16 = dB.carrFrac(i16);

% =====================================================================
% Figure 1  mechanism schematic
% =====================================================================
f = newfig(17, 7);
ax = axes('Position',[.06 .12 .88 .80]); hold(ax,'on'); axis(ax,[0 10 0 4]); axis off
rectangle('Position',[0 0 6 4],'FaceColor',[.94 .90 .86],'EdgeColor','none');
rectangle('Position',[6 0 4 4],'FaceColor',[.88 .94 .97],'EdgeColor','none');
plot([6 6],[0 4],'k--','LineWidth',1.2);
text(3,3.72,'host 5xFAD cortex','HorizontalAlignment','center','FontWeight','bold');
text(8,3.72,'wild-type graft','HorizontalAlignment','center','FontWeight','bold');
text(6,-0.22,'interface','HorizontalAlignment','center','FontSize',8);
rng(3);
for k = 1:16
    xc = 0.4+5.2*rand; yc = 0.45+3.0*rand; rr = 0.10+0.14*rand;
    rectangle('Position',[xc-rr yc-rr 2*rr 2*rr],'Curvature',[1 1], ...
              'FaceColor',[.80 .35 .30],'EdgeColor','none');
end
xx = linspace(6,10,200); lm = 0.9;
plot(xx, 1.05 + 0.95*exp(-(xx-6)/lm),'-','Color',[.20 .45 .80],'LineWidth',2);
text(7.0,2.45,'free A\beta, screened over {\itl}_m','Color',[.20 .45 .80],'FontSize',8.5);
annotation('arrow',[.52 .60],[.30 .30],'Color',[.20 .45 .80],'LineWidth',1.6);
for k = 1:3
    xc = 5.1 + 1.5*(k-1); yc = 0.75;
    rectangle('Position',[xc-.17 yc-.17 .34 .34],'Curvature',[1 1], ...
              'FaceColor',[.25 .55 .35],'EdgeColor','none');
    rectangle('Position',[xc-.06 yc-.06 .12 .12],'Curvature',[1 1], ...
              'FaceColor',[.80 .35 .30],'EdgeColor','none');
end
annotation('arrow',[.515 .645],[.235 .235],'Color',[.25 .55 .35],'LineWidth',1.8);
text(7.3,0.40,'carrier delivers a seed','Color',[.25 .55 .35],'FontSize',8.5);
text(0.25,0.35,{'microglia load at host plaques'},'FontSize',8.5);
annotation('textbox',[.655 .60 .33 .30],'String', ...
  {'\bf Serial mechanism','1.  carriers deliver the seed', ...
   '2.  diffusion delivers the mass'},'FitBoxToText','on', ...
   'BackgroundColor','w','EdgeColor',[.6 .6 .6],'FontSize',8.5);
savefigs(f, outdir, 'fig1');

% =====================================================================
% Figure 2  the seeding bottleneck
% =====================================================================
f = newfig(17, 6.5);
subplot(1,2,1); hold on
mm = logspace(-3,0,100);
S0 = p.detRatio * p.m_seed * p.c_h;   % surface of a just-visible deposit field
plot(mm, p.k_1*mm.^2, '-', 'Color',[.65 .25 .65]);
relRate = (1-p.rho)*p.k_r*mean(solB.b(i4,G.inGraft)) * ones(size(mm));
plot(mm, relRate, '-', 'Color',[.25 .55 .35]);
plot(mm, (p.k_a*mm + p.k_2*mm.^2)*S0, '-', 'Color',[.20 .45 .80]);
need = p.detRatio * p.m_seed * p.c_h / 28;
plot(mm, need*ones(size(mm)), 'k--');
set(gca,'XScale','log','YScale','log'); xlim([1e-3 1]); ylim([1e-12 1e2]);
xlabel('free A\beta, m (clean-tissue units)'); ylabel('initiation rate (units day^{-1})');
legend({'homogeneous nucleation k_1m^2','carrier release (1-\rho)k_rb', ...
        'growth on a detectable seed','rate needed by 4 weeks'}, ...
        'Location','southeast','Box','off','FontSize',7.5);
title('(a) routes to the first deposit');

subplot(1,2,2); hold on
jm = find(G.inGraft, 1, 'first') + round(0.5*sum(G.inGraft));
plot(solB.t/7, max(dB.areaPct(:,jm),1e-12), '-', 'Color',[.25 .55 .35]);
plot(solN.t/7, max(dN.areaPct(:,jm),1e-12), '-', 'Color',[.70 .40 .15]);
plot(xlim, [p.aThr p.aThr], 'k--');
set(gca,'YScale','log'); ylim([1e-9 1e1]); xlabel('weeks after transplantation');
ylabel('deposit at mid-graft depth (% area)');
legend({'carriers active','carriers inert','detection threshold'}, ...
       'Location','southeast','Box','off','FontSize',7.5);
title('(b) no carriers, no detectable deposit');
savefigs(f, outdir, 'fig2');

% =====================================================================
% Figure 3  space-time maps
% =====================================================================
f = newfig(17, 6);
names = {'microglia  \nu','carrier cargo  b','deposited A\beta  F'};
flds  = {solB.nu, solB.b, solB.F};
for k = 1:3
    subplot(1,3,k);
    Z = flds{k};
    imagesc(G.rc, solB.t/7, Z); set(gca,'YDir','normal'); hold on
    plot([p.Rg p.Rg], [0 max(solB.t)/7], 'w--', 'LineWidth',1.2);
    tt = linspace(0, max(solB.t), 50);
    plot(max(p.Rg - dB.Vmig*tt, 0), tt/7, 'w-', 'LineWidth',1.0);
    plot(max(p.Rg - R.ell_m_host, 0)*[1 1], [0 max(solB.t)/7], 'w:', 'LineWidth',1.0);
    xlim([0 1.2]); xlabel('r (mm)'); if k==1, ylabel('weeks'); end
    title(sprintf('(%c) %s', 'a'+k-1, names{k})); colorbar
end
colormap(parula);
savefigs(f, outdir, 'fig3');

% =====================================================================
% Figure 4  calibration
% =====================================================================
fprintf('running lesion geometry ...\n');
Lw = p210_lesion('baseline');  Lk = p210_lesion('irf8ko');
R.lesNuRatio = Lw.nuWin(end)/Lw.nuWin(1);
R.lesAbWT    = Lw.AbWinRel(end);
R.lesAbKO    = Lk.AbWinRel(end);
% The genotype comparison must be made on BASELINE-NORMALIZED increments.
% The 15.0-fold ratio of absolute 24 h areas is dominated by a difference
% already present at t = 0 (16.35 vs 1.58 um^2, implied by Fig. 3c with
% Fig. 3d); the model starts both genotypes level, so only the increment
% ratio tests it.  See D.lesAb.t0Ratio and D.lesAb.incRatio.
R.lesAbRatio = Lw.AbWinRel(end) / max(Lk.AbWinRel(end), eps);

f = newfig(17, 6.5);
subplot(1,3,1); hold on
plot(Lw.tHr, Lw.nuWin/Lw.nuWin(1), '-', 'Color',[.25 .55 .35]);
plot(Lk.tHr, Lk.nuWin/Lk.nuWin(1), '-', 'Color',[.70 .40 .15]);
errorbar(D.lesNu.tHr, D.lesNu.wt/D.lesNu.wt(1), D.lesNu.wtErr/D.lesNu.wt(1), ...
         'ko','MarkerFaceColor','k','MarkerSize',4,'LineStyle','none');
plot(D.lesNu.tHr, D.lesNu.ko/D.lesNu.ko(1), 'ks','MarkerSize',4,'LineStyle','none');
xlabel('hours after lesion'); ylabel('microglia in window (relative)');
ylim([0.5 4.2]);          % headroom so the legend clears the 24 h error bar
legend({'Irf8^{+/+} model','Irf8^{-/-} model','Fig. 3e Irf8^{+/+}','Fig. 3e Irf8^{-/-}'}, ...
       'Location','northwest','Box','off','FontSize',7);
title('(a) carrier recruitment  [FIT]','FontSize',8.5);

subplot(1,3,2); hold on
plot(Lw.tHr, Lw.AbWinRel, '-', 'Color',[.25 .55 .35]);
plot(Lk.tHr, Lk.AbWinRel, '-', 'Color',[.70 .40 .15]);
xlabel('hours after lesion'); ylabel('A\beta in window (relative)');
ylim([1 6.2]);            % headroom for the legend
legend({'Irf8^{+/+}','Irf8^{-/-}'},'Location','northwest','Box','off','FontSize',7.5);
title('(b) cargo arrival  [PREDICTION]','FontSize',8.5);

subplot(1,3,3); hold on
mi = [Lw.AbWinRel(end) Lk.AbWinRel(end)];
di = [D.lesInc.wt(3)   D.lesInc.ko(3)];
ei = [D.lesInc.wtErr(3) D.lesInc.koErr(3)];
bar(1:2, mi, 0.6, 'FaceColor',[.75 .80 .85]);
errorbar(1:2, di, ei, 'ko','MarkerFaceColor','k','MarkerSize',4,'LineStyle','none');
set(gca,'XTick',1:2,'XTickLabel',{'Irf8^{+/+}','Irf8^{-/-}'});
ylabel('A\beta in window at 24 h (x baseline)');
ylim([0 20]);
title('(c) 24-h increment  [PREDICTION]','FontSize',8.5);
text(0.55, 0.94, sprintf('ratio %.2f (model)', R.lesAbRatio), 'Units','normalized','FontSize',7.5);
text(0.55, 0.87, sprintf('%.2f (measured)', D.lesAb.incRatio), 'Units','normalized','FontSize',7.5);
savefigs(f, outdir, 'fig4');

% =====================================================================
% Figure 5  validation of the grafting experiments
% =====================================================================
f = newfig(17, 11);
subplot(2,2,1); hold on
plot(dB.t/7, dB.nuGraft, '-', 'Color',[.25 .55 .35]);
errorbar(D.graftNu.t/7, D.graftNu.val, D.graftNu.err, 'ko', ...
         'MarkerFaceColor','k','MarkerSize',4,'LineStyle','none');
xlabel('weeks'); ylabel('graft microglia (mm^{-2})');
legend({'model','Fig. 1d'},'Box','off','Location','northwest','FontSize',7.5);
title('(a) graft infiltration  [FITTED]');

subplot(2,2,2); hold on
plot(dB.t/7, dB.loadPct, '-', 'Color',[.25 .55 .35]);
plot(28/7, D.irf8Load.val(1), 'ks','MarkerFaceColor',[.9 .9 .2],'MarkerSize',7);
errorbar(D.graftTime.t/7, D.graftTime.val, D.graftTime.err, 'ko', ...
         'MarkerFaceColor','w','MarkerSize',4,'LineStyle','none');
xlabel('weeks'); ylabel('graft A\beta load (% area)');
legend({'model','Fig. 2h (calibration)','ED Fig. 1c (other cohort)'}, ...
       'Box','off','Location','northwest','FontSize',7.5);
title('(b) graft amyloid time course');

% Four perturbations, three of which fit no amyloid parameter.  The measured
% bars for BLZ945, aged and Cx3cr1-/- are their published within-cohort
% ratios applied to the Fig. 2h control, because each comes from a cohort
% with its own absolute level.
dC  = p210_post(p210_solve(p210_params('cx3cr1ko'), TS));
lab = {'Irf8^{+/-}','Irf8^{-/-}','BLZ945','aged','Cx3cr1^{-/-}'};
nb  = numel(lab);
subplot(2,2,3); hold on
mv = [dB.loadPct(i4) dK.loadPct(i4) dD.loadPct(i4) dA.loadPct(i4) dC.loadPct(i4)];
bar(1:nb, mv, 0.6, 'FaceColor',[.75 .80 .85]);
dv = [D.irf8Load.val(1) D.irf8Load.val(2) ...
      D.irf8Load.val(1)*D.deplete.rel D.irf8Load.val(1)*D.aged.rel ...
      D.irf8Load.val(1)*D.cx3cr1.rel];
de = [D.irf8Load.err(1) D.irf8Load.err(2) ...
      D.irf8Load.val(1)*D.deplete.relErr D.irf8Load.val(1)*D.aged.relErr ...
      D.irf8Load.val(1)*D.cx3cr1.relErr];
errorbar(1:nb, dv, de, 'ko','MarkerFaceColor','k','MarkerSize',4,'LineStyle','none');
set(gca,'XTick',1:nb,'XTickLabel',lab); xtickangle(20);
ylabel('graft A\beta at 4 weeks (% area)');
title('(c) perturbations at 4 weeks');

subplot(2,2,4); hold on
bar(1:nb, mv./mv(1), 0.6, 'FaceColor',[.75 .80 .85]);
errorbar(1:nb, dv./dv(1), de./dv(1), 'ko','MarkerFaceColor','k','MarkerSize',4,'LineStyle','none');
set(gca,'XTick',1:nb,'XTickLabel',lab); xtickangle(20);
ylabel('relative to Irf8^{+/-}'); ylim([0 1.3]);
title('(d) the same, normalised');
savefigs(f, outdir, 'fig5');

% =====================================================================
% Figure 6  the two signatures and the apportionment
% =====================================================================
fprintf('running depletion sweep ...\n');
fr  = linspace(0, 0.95, 9);        % fraction of microglia removed
Ldep = zeros(size(fr)); Hdep = zeros(size(fr));
for k = 1:numel(fr)
    pk = p210_params('baseline'); pk.N = 240; pk.nu0 = max(1-fr(k), 1e-3);
    dk = p210_post(p210_solve(pk, [0 7 14 21 28]));
    Ldep(k) = dk.loadPct(end); Hdep(k) = dk.hostPct(end);
end
R.deplCurve = [fr(:) Ldep(:) Hdep(:)];
R.depl80 = interp1(fr, Ldep, 0.80) / Ldep(1);

f = newfig(17, 6.5);
subplot(1,3,1); hold on
cols = lines(4); tt = [14 28 56 112];
for k = 1:4
    j = find(TS==tt(k));
    plot(p.Rg - G.rc(G.inGraft), dB.areaPct(j,G.inGraft), '-', 'Color',cols(k,:));
end
dd = linspace(0, p.Rg, 50);
j = find(TS==112);
F0 = dB.areaPct(j, find(G.inGraft,1,'last'));
plot(dd, F0*exp(-dd/R.ell_m_graft), 'k--');
xlabel('depth into graft (mm)'); ylabel('deposit (% area)');
legend({'2 wk','4 wk','8 wk','16 wk','exp(-x/{\itl}_m)'},'Box','off','FontSize',7.5);
title('(a) amplitude is diffusion-set');

subplot(1,3,2); hold on
ok = ~isnan(dB.onsetT);
plot(dB.depth(ok), dB.onsetT(ok)/7, 'o-', 'Color',[.25 .55 .35],'MarkerSize',3);
ok2 = ~isnan(dK.onsetT);
if any(ok2), plot(dK.depth(ok2), dK.onsetT(ok2)/7, 'o-','Color',[.70 .40 .15],'MarkerSize',3); end
xlabel('depth into graft (mm)'); ylabel('onset of detectable deposit (weeks)');
legend({'Irf8^{+/-}','Irf8^{-/-}'},'Box','off','Location','northwest','FontSize',7.5);
title('(b) timing is carrier-set');

subplot(1,3,3); hold on
plot(fr, Ldep/Ldep(1), 'o-', 'Color',[.25 .55 .35],'MarkerSize',3);
errorbar(1 - D.deplete.nuRel, D.deplete.rel, D.deplete.relErr, 'ks', ...
         'MarkerFaceColor',[.9 .9 .2],'MarkerSize',7);
xlabel('fraction of microglia removed'); ylabel('graft A\beta, relative to intact');
legend({'model','BLZ945, Ext. Data Fig. 7'},'Box','off','Location','southwest','FontSize',7.5);
title('(c) depletion  [BLZ945 FITTED]');
savefigs(f, outdir, 'fig6');

% =====================================================================
% Figure 7  clearance against dissemination
% =====================================================================
fprintf('running export-length and graft-size sweeps ...\n');
krv = logspace(-3, 0, 9);
Lam = zeros(size(krv)); dGraft = zeros(size(krv)); dHost = zeros(size(krv));
pn = p210_params('nocarrier'); pn.N = 240;
dn0 = p210_post(p210_solve(pn, [0 28]));
for k = 1:numel(krv)
    pk = p210_params('baseline'); pk.N = 240; pk.k_r = krv(k);
    dk = p210_post(p210_solve(pk, [0 28]));
    Lam(k)    = dk.Lambda;
    dGraft(k) = dk.loadPct(end) - dn0.loadPct(end);
    dHost(k)  = dk.hostPct(end) - dn0.hostPct(end);
end
R.LamSweep = [Lam(:) dGraft(:) dHost(:)];

Rgv = [0.15 0.25 0.35 0.45 0.60 0.80 1.10];
Lrg = zeros(size(Rgv));
for k = 1:numel(Rgv)
    pk = p210_params('baseline'); pk.N = 240; pk.Rg = Rgv(k);
    dk = p210_post(p210_solve(pk, [0 28]));
    Lrg(k) = dk.loadPct(end);
end
R.RgSweep = [Rgv(:) Lrg(:)];

f = newfig(17, 6.5);
subplot(1,3,1); hold on
plot(Lam, dGraft, 'o-','Color',[.25 .55 .35],'MarkerSize',3);
plot(Lam, dHost,  'o-','Color',[.70 .40 .15],'MarkerSize',3);
plot(xlim,[0 0],'k:');
set(gca,'XScale','log'); xlabel('\Lambda = {\itl}_x / {\itl}_m');
ylabel('change in A\beta due to carriers (% area)');
legend({'graft (naive tissue)','host (plaque-bearing)'},'Box','off','Location','northwest','FontSize',7.5);
title('(a) net effect of carriers');

subplot(1,3,2); hold on
plot(fr, Ldep-Ldep(1), 'o-','Color',[.25 .55 .35],'MarkerSize',3);
plot(fr, Hdep-Hdep(1), 'o-','Color',[.70 .40 .15],'MarkerSize',3);
plot(xlim,[0 0],'k:');
xlabel('fraction of microglia removed'); ylabel('change in A\beta (% area)');
legend({'graft','host'},'Box','off','Location','northwest','FontSize',7.5);
title('(b) depletion: opposite signs');

subplot(1,3,3); hold on
plot(Rgv, Lrg, 'o-','Color',[.25 .55 .35],'MarkerSize',3);
rr = linspace(min(Rgv),max(Rgv),50);
plot(rr, Lrg(4)*(Rgv(4)./rr), 'k--');
xlabel('graft radius R_g (mm)'); ylabel('whole-graft A\beta at 4 weeks (% area)');
legend({'model','{\itl}_m / R_g scaling'},'Box','off','FontSize',7.5);
title('(c) graft-size confound');
savefigs(f, outdir, 'fig7');

% =====================================================================
% Figure 8  operating diagram
% =====================================================================
fprintf('running the operating diagram ...\n');
% The two axes are the two independent things the experiments do to the
% carrier population: how many carriers reach the graft, and how much each
% one picks up.  Irf8-/- and BLZ945 move the first at fixed uptake; aged
% recipients move the second at an unchanged density.  The two published
% axes of perturbation are therefore close to orthogonal in this plane.
nk = 7;
phg = logspace(-1, 0.15, nk);
nug = logspace(-1.2, 0.1, nk);
Zl = zeros(nk); Zc = zeros(nk);
p0 = p210_params('baseline');
Dd = p210_data();
for a = 1:nk
  for b2 = 1:nk
    pk = p210_params('baseline'); pk.N = 200;
    pk.phi = phg(a);  pk.k_u = p0.k_u * phg(a);
    pk.nu0 = nug(b2);
    dk = p210_post(p210_solve(pk, [0 28]));
    Zl(b2,a) = dk.loadPct(end);
    Zc(b2,a) = min(dk.carrFrac(end), 1);
  end
end
R.opLoad = Zl;  R.opCarr = Zc;  R.opPhi = phg;  R.opNu = nug;

f = newfig(17, 6.5);
subplot(1,2,1);
contourf(log10(phg), log10(nug), Zl, 12, 'LineColor','none'); hold on
plot(0,                      0,                       'wp','MarkerFaceColor','w','MarkerSize',11);
plot(0,                      log10(Dd.irf8Nu.rel),    'wo','MarkerFaceColor','none','MarkerSize',8,'LineWidth',1.4);
plot(0,                      log10(Dd.deplete.nuRel), 'ws','MarkerFaceColor','none','MarkerSize',8,'LineWidth',1.4);
plot(log10(p0.phiAged),      log10(Dd.aged.nuRel),    'wd','MarkerFaceColor','none','MarkerSize',8,'LineWidth',1.4);
xlabel('log_{10} (uptake per carrier / control)'); ylabel('log_{10} (\nu_0 / \nu_0^{ref})');
title('(a) graft A\beta at 4 weeks (%)'); colorbar
text(0,                 0.06,                        ' control',   'Color','w','FontSize',8);
text(0,                 log10(Dd.irf8Nu.rel)+0.06,   ' Irf8^{-/-}','Color','w','FontSize',8);
text(0,                 log10(Dd.deplete.nuRel)-0.13,' BLZ945',    'Color','w','FontSize',8);
text(log10(p0.phiAged), log10(Dd.aged.nuRel)+0.06,   ' aged',      'Color','w','FontSize',8);

subplot(1,2,2);
contourf(log10(phg), log10(nug), Zc, 12, 'LineColor','none'); hold on
plot(0, 0, 'wp','MarkerFaceColor','w','MarkerSize',11);
xlabel('log_{10} (uptake per carrier / control)'); ylabel('log_{10} (\nu_0 / \nu_0^{ref})');
title('(b) carrier-released fraction'); colorbar
savefigs(f, outdir, 'fig8');

% =====================================================================
% Supplementary S1  field profiles
% =====================================================================
f = newfig(17, 11);
lab = {'free A\beta  m','deposit  F (% area)','cargo  b','microglia  \nu'};
dat = {solB.m, dB.areaPct, solB.b, solB.nu};
for k = 1:4
    subplot(2,2,k); hold on
    for j = [i2 i4 i16]
        plot(G.rc, dat{k}(j,:), '-');
    end
    plot([p.Rg p.Rg], ylim, 'k--');
    xlim([0 1.5]); xlabel('r (mm)'); ylabel(lab{k});
    if k==1, legend({'2 wk','4 wk','16 wk','interface'},'Box','off','FontSize',7.5); end
end
savefigs(f, outdir, 'figS1');

% =====================================================================
% Supplementary S2  numerical verification
% =====================================================================
fprintf('running grid convergence ...\n');
Ns = [160 320 640];  conv = zeros(numel(Ns),2);
for k = 1:numel(Ns)
    pk = p210_params('baseline'); pk.N = Ns(k);
    dk = p210_post(p210_solve(pk, [0 28 112]));
    conv(k,:) = [dk.loadPct(2) dk.loadPct(3)];
end
R.conv = [Ns(:) conv];
f = newfig(17, 6.5);
subplot(1,2,1); hold on
plot(Ns, conv(:,1)/conv(end,1), 'o-'); plot(Ns, conv(:,2)/conv(end,2), 's-');
xlabel('grid points N'); ylabel('graft load, relative to N = 640');
legend({'4 weeks','16 weeks'},'Box','off','FontSize',7.5);
title('(a) grid convergence');
subplot(1,2,2); hold on
plot(G.rc, min(solB.nu,[],1), '-'); plot(G.rc, min(solB.b,[],1), '-');
plot(xlim,[0 0],'k:'); xlim([0 1.5]);
xlabel('r (mm)'); ylabel('minimum over time');
legend({'\nu','b'},'Box','off','FontSize',7.5);
title('(b) positivity of the upwind scheme');
savefigs(f, outdir, 'figS2');

% =====================================================================
% Supplementary S3  one-at-a-time sensitivity
% =====================================================================
fprintf('running sensitivity ...\n');
names = {'D_n','chi_s','k_p','k_r','k_d','b_max','D_m','lam_m','k_a','k_2','Rg','g_n'};
sens = zeros(numel(names),2);
for k = 1:numel(names)
    for s = 1:2
        pk = p210_params('baseline'); pk.N = 200;
        f0 = 0.5*(s==1) + 2.0*(s==2);
        pk.(names{k}) = pk.(names{k}) * f0;
        dk = p210_post(p210_solve(pk, [0 28]));
        sens(k,s) = dk.loadPct(end);
    end
end
R.sens = sens;  R.sensNames = names;
f = newfig(12, 11);
base4 = R.load4;
barh(1:numel(names), [sens(:,1)-base4, sens(:,2)-base4], 'stacked');
set(gca,'YTick',1:numel(names),'YTickLabel',strrep(names,'_','\_'));
xlabel('change in graft A\beta at 4 weeks (% area)');
legend({'halved','doubled'},'Box','off','Location','southeast','FontSize',7.5);
title('one-at-a-time sensitivity');
savefigs(f, outdir, 'figS3');

fprintf('\nall figures written to %s\n', outdir);
end

% =========================================================================
function f = newfig(wcm, hcm)
f = figure('Units','centimeters','Position',[2 2 wcm hcm], ...
           'Color','w','PaperPositionMode','auto');
end

function savefigs(f, outdir, name)
print(f, fullfile(outdir,[name '.png']), '-dpng', '-r300');
try, savefig(f, fullfile(outdir,[name '.fig'])); catch, end
close(f);
end
