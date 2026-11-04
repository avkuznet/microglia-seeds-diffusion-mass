function D = p210_data()
%P210_DATA  Measurements of d'Errico et al. (2022), from the published
%   Source Data spreadsheets that accompany the paper.
%
%   These are the authors' own per-animal values (Source Data Fig. 1-3 and
%   Extended Data Fig. 4, 7), not values digitised from printed bar graphs.
%   Means and s.e.m. are as computed in those files.  Each entry records the
%   figure panel it comes from.
%
%   NOTE ON THE CONTROL GRAFT LOAD.  Three cohorts give a 4-week control
%   graft load: Fig. 2h (Irf8+/-)  0.586 +/- 0.093 % area, n = 14
%              ED Fig. 4j (young)  1.083 +/- 0.231 % area, n = 12
%              ED Fig. 7g (vehicle) 0.274 +/- 0.067 % area, n = 13
%   They span a factor of four.  The calibration uses the Fig. 2h cohort,
%   which is the one the Irf8 comparison is made within; the spread between
%   cohorts is the dominant uncertainty on any absolute load.

% Fig. 1d  microglial density inside grafts, cells/mm^2 on 25 um sections
D.graftNu.t   = [14 28];
D.graftNu.val = [750.68 516.81];
D.graftNu.err = [ 83.87  31.55];
D.graftNu.n   = [7 13];

% Extended Data Fig. 1c  graft amyloid load against time, one cohort.
%   The only published 16-week graft load.  Note the cohort offset: this
%   cohort's 4-week value is 1.083, not the 0.586 of Fig. 2h that the model
%   is calibrated to, so the within-cohort ratios are the meaningful
%   comparison and the absolute values are not.
D.graftTime.t    = [14 28 112];           % days
D.graftTime.val  = [0.1975 1.083 1.924];  % % of area
D.graftTime.err  = [0.0754 0.2312 0.5186];
D.graftTime.n    = [7 12 12];            % GRAFTS
D.graftTime.nMice = [4 9 8];             % MICE, from the ED Fig. 1c caption
%   Every graft-level panel in d'Errico et al. is captioned "each symbol
%   represents one graft from n mice", so all of these counts, and all the
%   standard errors computed from them, are clustered: Fig. 2h is 14 grafts
%   from 7 mice and 12 from 6, ED Fig. 5c is 12 from 9 and 10 from 7, and
%   ED Fig. 7g is 13 from 8 and 17 from 8.  The published ED Fig. 1c test,
%   one-way ANOVA with F(2,28), has 31 total observations, so the source
%   analysis is itself at the graft level.
%   The individual animals, needed because these distributions are strongly
%   right-skewed and a ratio of means with Gaussian errors misrepresents
%   them; see p210_ratiostats.  Two of the seven 2-week grafts carry no
%   detectable amyloid at all, which is why the 2-week point admits no
%   logarithmic treatment.
D.graftTime.v2   = [0 0.262454943 0 0.382880457 0.091065120 0.522132039 0.124150485];
D.graftTime.v4   = [0.546905206 0.641315279 1.557610178 0.360934505 1.871460231 ...
                    0.274481652 1.431207773 1.717974436 0.605357143 0.036588748 ...
                    2.719512385 1.238226012];
D.graftTime.v16  = [1.423117836 6.663614294 0.361349238 2.392914338 2.584404519 ...
                    0.961868966 3.729559147 0.137505262 0.838853247 1.103975133 ...
                    1.454696200 1.436451307];

% Fig. 2g  graft microglial density, Irf8+/- against Irf8-/-
D.irf8Nu.val  = [519.50 68.89];
D.irf8Nu.err  = [ 34.16 21.52];
D.irf8Nu.rel  = 68.89 / 519.50;        % 0.133

% Fig. 2h  graft amyloid load, % of area, Irf8+/- against Irf8-/-
D.irf8Load.val = [0.5863 0.1433];
D.irf8Load.err = [0.0932 0.0276];
D.irf8Load.rel = 0.1433 / 0.5863;      % 0.244
D.irf8Load.relErr = D.irf8Load.rel * hypot(0.0276/0.1433, 0.0932/0.5863);

% Fig. 3e  microglia counted in the 150 x 150 um lesion window
D.lesNu.tHr   = [0 6 24];
D.lesNu.wt    = [13.33 15.33 24.25];
D.lesNu.wtErr = [ 2.58  4.32 11.67];
D.lesNu.ko    = [11.00  8.20 13.00];
D.lesNu.rel24 = 24.25 / 13.33;         % 1.82

% Fig. 3c  amyloid area at the lesion at 24 h, um^2, in the same window
D.lesAb.val   = [153.70 10.22];
D.lesAb.err   = [ 45.62  1.386];
D.lesAb.ratio = 153.70 / 10.22;        % 15.0
D.lesAb.winUm2 = 150*150;              % window area, um^2

% Fig. 3d  relative increment of amyloid at the lesion site
D.lesInc.tHr = [0 6 24];
D.lesInc.wt  = [1 2.938 9.402];
D.lesInc.ko  = [1 1.020 6.452];
D.lesInc.wtErr = [0 1.444 8.588];
D.lesInc.koErr = [0 0.249 6.235];
% Initial amyloid in the window implied by Fig. 3c and Fig. 3d, % of area:
D.lesWin0Pct = (153.70/9.402) / D.lesAb.winUm2 * 100;   % 0.0727
%   DECOMPOSITION OF THE 15-FOLD GENOTYPE DIFFERENCE AT 24 h.  Fig. 3c gives
%   absolute areas at 24 h and Fig. 3d the increment each genotype made over
%   its own starting level, so the starting levels can be recovered and they
%   are not the same: 16.35 um^2 in Irf8+/+ against 1.58 in Irf8-/-.  The
%   15.0-fold ratio of the 24 h areas is therefore 10.3 already present
%   before any delivery, times 1.46 contributed by delivery.  A model that
%   starts both genotypes at the same window level can only be tested
%   against the second factor.
D.lesAb.t0   = D.lesAb.val ./ [D.lesInc.wt(3) D.lesInc.ko(3)];   % 16.35, 1.58
D.lesAb.t0Ratio  = D.lesAb.t0(1) / D.lesAb.t0(2);                % 10.3
D.lesAb.incRatio = D.lesInc.wt(3) / D.lesInc.ko(3);              % 1.46

% Extended Data Fig. 7f,g  BLZ945 depletion: graft microglia and graft load
D.deplete.nuRel = 77.37 / 585.90;      % 0.132 (measured, not the stated 80%)
D.deplete.rel   = 0.06994 / 0.27350;   % 0.256
D.deplete.val   = [0.27350 0.06994];
D.deplete.err   = [0.06719 0.01493];
D.deplete.relErr = D.deplete.rel * hypot(0.01493/0.06994, 0.06719/0.27350);

% Extended Data Fig. 4i,j  aged recipients
D.aged.nuRel = 556.85 / 516.81;        % 1.08: densities equivalent (P = 0.54)
D.aged.rel   = 0.192 / 1.083;          % 0.177
D.aged.val   = [1.083 0.192];
D.aged.err   = [0.2312 0.1391];
D.aged.relErr = D.aged.rel * hypot(0.1391/0.192, 0.2312/1.083);

% Extended Data Fig. 6b,d  HOST cortical and hippocampal plaque density.
%   This is the deposit NUMBER density, and it is the quantity the discrete
%   observation operator needs.  Note the second column against the third:
%   removing Irf8, and with it seven microglia in eight, leaves the host
%   deposit density unchanged.  Only naive tissue depends on the carriers.
D.plaqDens.ctx     = [54.05 38.57 41.18];     % +/+, +/-, -/-  per mm^2
D.plaqDens.ctxErr  = [15.07  5.67  5.911];
D.plaqDens.hipp    = [23.14 16.08 13.69];
D.plaqDens.rel     = 41.18 / 38.57;           % 1.068: unchanged
D.plaqDens.sectMm  = 0.025;                   % section thickness
D.plaqDens.diaMm   = 0.020;                   % representative plaque size
% Abercrombie: per mm^3 = per mm^2 / (section + object diameter)
D.plaqDens.perMm3  = 38.57 / (0.025 + 0.020); % 857 per mm^3

% Extended Data Fig. 6f  microglia associated with each cortical plaque.
%   Used as the cortical microglial population of an Irf8-/- recipient
%   relative to control, the deposit density being the same in both.
D.irf8Cortex.val = [2.744 2.504 0.890];       % +/+, +/-, -/-
D.irf8Cortex.err = [0.2544 0.09543 0.09397];
D.irf8Cortex.rel = 0.890 / 2.504;             % 0.355

% Extended Data Fig. 6h  amyloid carried per microglial cell, cortex
D.irf8Cargo.val = [7.174 4.753];              % Irf8+/-, Irf8-/-
D.irf8Cargo.err = [1.064 0.7294];
D.irf8Cargo.rel = 4.753 / 7.174;              % 0.663

% Fig. 2d,e  process motility of Irf8-/- microglia, 0-30 um from the lesion
D.irf8Motility.dist  = [15.615 10.642];       % um travelled, +/+ vs -/-
D.irf8Motility.speed = [0.446  0.304];        % um/min
D.irf8Motility.rel   = mean([10.642/15.615, 0.304/0.446]);   % 0.681

% Extended Data Fig. 5b,c  Cx3cr1-/- recipients: the same number of
%   microglia, a quarter of the amyloid.  A per-cell effect, not a number
%   effect, and so a test of a different axis of the model from Irf8/BLZ945.
D.cx3cr1.nuRel = 506.1 / 516.8;               % 0.979
D.cx3cr1.val   = [1.083 0.2961];              % % area, control then Cx3cr1-/-
D.cx3cr1.err   = [0.2312 0.0772];
D.cx3cr1.rel   = 0.2961 / 1.083;              % 0.273
D.cx3cr1.relErr = D.cx3cr1.rel * hypot(0.0772/0.2961, 0.2312/1.083);

% Extended Data Fig. 3b,c  Ccr2-/- recipients (peripheral monocytes)
D.ccr2.nuRel = 511.1 / 707.4;                 % 0.722
D.ccr2.rel   = 0.4531 / 0.2072;               % 2.19, with wide scatter

% ---- load against carrier number -----------------------------------------
%   Two independent perturbations of the carrier number, one genetic and one
%   pharmacological, give the same power law between the graft load and the
%   microglial density:  load ~ nu^gamma.
D.gamma.val = 0.687;                          % weighted mean of the two
D.gamma.err = 0.123;
D.gamma.each = [log(0.2444)/log(0.1326), log(0.2557)/log(0.1321)];  % .70 .67

% Extended Data Fig. 4c,d  cargo carried per aged microglial cell
D.agedCargo.pctCells  = 30.533 / 51.285;   % 0.595  A-beta-associated cells
D.agedCargo.perCell   = 0.5139 / 1.0410;   % 0.494  A-beta clusters per cell
D.agedCargo.rel       = 0.494;             % used as the measured uptake factor
end
