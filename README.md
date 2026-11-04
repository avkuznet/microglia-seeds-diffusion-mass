# P210 — microglial carrier–diffusion model of Aβ propagation into naive tissue

MATLAB implementation for *"Microglia deliver the seeds and diffusion delivers the mass:
a transport model of amyloid-β propagation into naive brain tissue."*

Base MATLAB only, no toolboxes; developed and run on R2023a. MIT licence.

## Quick start

```matlab
p210_calibrate      % fit the six constants, writes p210_fitted.mat   (~4 min)
run_all             % print every number quoted in the paper          (~30 s)
make_figures        % write Figs 1-8 and S1-S3 to ../figures           (~60 s)
```

`run_all('quick')` skips calibration and uses the stored fit in `p210_fitted.mat`.
The console log saved as `run_all_output.txt` is from a quick run of that kind; its
header records which it was.

## The model

Six fields on a spherically symmetric domain, with the graft occupying `0 ≤ r < Rg`
and the host `Rg < r ≤ Rout`:

| field | meaning |
|---|---|
| `m`  | free, diffusible Aβ |
| `F`  | deposited, immobile fibrillar Aβ |
| `b`  | microglia-internalised Aβ (cargo) |
| `nu` | microglial density, in multiples of homeostatic |
| `s`  | injury-derived chemoattractant |
| `c`  | number density of distinguishable deposits, per mm³ |

All three amyloid pools share one unit — the free concentration in clean tissue — so
mass is conserved between them: what a carrier takes up is removed from the pool it
came from, and what it releases is returned. Cargo is advected by the *same* discrete
operator as the cells, which is what keeps the per-carrier load `b/nu` free of
spurious extrema at the migration front.

## Files

| file | purpose |
|---|---|
| `p210_params.m` | parameter sets: `'baseline'`, `'irf8ko'`, `'depleted'`, `'aged'`, `'nocarrier'` |
| `p210_grid.m` | finite-volume grid (planar, cylindrical or spherical) and Jacobian sparsity |
| `p210_rhs.m` | right-hand side: diffusion, chemotaxis, uptake, release, deposition |
| `p210_solve.m` | method-of-lines driver around `ode15s` |
| `p210_post.m` | observation operator and diagnostics: loads, densities, onset times, mass budget, length scales |
| `p210_scales.m` | transport length scales `ell_m`, `ell_x` and `Lambda` on their own fine time grid; the peak chemotactic drift occurs in the first few days and is understated by a factor of two if read off the published output times |
| `p210_ratiostats.m` | compares a model ratio with a measured one three ways: ratio of means with Gaussian errors, log scale, and bootstrap. Graft loads are strongly right-skewed, so the first misrepresents them and the log scale is the one quoted |
| `p210_derive.m` | quantities that follow from the declared ones (`kap`, `F_h`, `c_h`, `m_c`) |
| `p210_lesion.m` | focal laser-lesion geometry used for calibration and the out-of-sample test |
| `p210_calibrate.m` | fits the six free constants to six published targets |
| `p210_data.m` | measurements from the Source Data files of d'Errico et al. (2022) |
| `run_all.m` | reproduces every number quoted in the paper |
| `make_figures.m` | Figs 1–8 and S1–S3 |
| `p210_competition.m` | the two observation operators recalibrated and compared |
| `p210_diag.m` | scratch diagnostic for the monomer balance |
| `smoke.m` | smoke test: solves the baseline case on a coarse time grid and prints loads, densities, carrier fraction and field minima |
| `p210_fitted.mat` | stored calibration |

## What is fitted and what is not

Six constants, one target each. Five are taken from the **Source Data** files published
with d'Errico et al. (2022) rather than digitised from the printed figures; the sixth,
`k_a`, is fitted to independent interstitial-Aβ depletion measurements.
Everything else is measured, derived from a cited measurement, or declared, and every
other comparison in the paper is a prediction.

| constant | value | fitted to |
|---|---|---|
| `k_a` | 0.2822 day⁻¹ | host interstitial Aβ depleted 80% at 4 weeks |
| `perPercent` | 34.49 units per 1% area | graft load 0.586% at 4 weeks — a unit conversion only |
| `chi_s` | 0.0835 mm² day⁻¹ per unit s | graft microglia 751 mm⁻² at 2 weeks (Fig. 1d) |
| `fnu_graft` | 0.274 | graft microglia 517 mm⁻² at 4 weeks (Fig. 1d) |
| `m_seed` | 4.62 × 10⁻⁴ units per deposit | BLZ945 load is 0.256 of vehicle (ED Fig. 7f,g) |
| `Q0les` | 107.5 day⁻¹ | microglia at a laser lesion rise 1.82-fold in 24 h (Fig. 3e) |

`chi_s` and `fnu_graft` fix the microglial field, which is an input to the amyloid
problem rather than a result of it.

**The aged arm fits nothing either.** Its uptake per cell is the measured 0.494 of
ED Fig. 4d and its host load is declared, so its graft load is a prediction. Earlier
versions also fitted a reduced cargo release rate for aged carriers, because with a
mass-proportional readout the measured two-fold uptake deficit fell far short of the
observed ratio — it needed a thirty-fold loss. With the discrete-deposit readout the
load follows delivery sublinearly and the measured deficit alone lands inside the
experimental error. `p.useAgedRelease` restores the old fit, which moves the predicted
ratio from 0.26 to 0.18 against a measured 0.177 ± 0.134; Njie et al. (2012) support
that mechanism independently (young microglia expel internalised Aβ within about three
hours without degrading it, aged microglia show no significant expulsion), but the data
no longer require it.

`m_seed` has an independent estimate it was not fitted to. A loaded microglial cell
carries about one Aβ cluster (ED Fig. 4d counts 1.04 per cell) and leaves it behind
when it releases, so one release event should start one deposit of mass
`b_max/n0_3D` = 9.4 × 10⁻⁴ units. The fitted value is 0.49 of that.

**The Irf8⁻ᐟ⁻ arm fits nothing.** Its graft carrier density is imposed at the measured
0.133 of control (Fig. 2g), exactly as BLZ945's is, so the amyloid load is a prediction.

Measured inputs: `D_m` = 62.3 µm² s⁻¹ (Waters 2010); `lam_m` = 9.8 × 10⁻⁵ s⁻¹,
half-life 1.97 h (Cirrito et al. 2003); host deposit number density `c_h` = 857 mm⁻³,
from the 38.6 ± 5.7 cortical plaques per mm² of ED Fig. 6b with an Abercrombie
correction. Deposit levels are declared as **areas**.

## The observation operator

Amyloid is laid down as discrete deposits and the histological readout is the area
they cover, not the mass they contain. The model carries their number density `c` as a
sixth field and reads

```
area  =  kap * c^(1/3) * F^(2/3) * vis(F/c)
```

— `c` deposits each of mass `F/c`, each showing an area proportional to its mass to the
two-thirds, each visible once it has grown past about ten times the seed it started
from. The mass-proportional operator of the first version is the special case in which
the deposit number is pinned at the host value; `p.readout = 'meanfield'` restores it,
and `p210_competition` recalibrates both and compares them.

This is where the reviewers' request for sub-voxel competition belongs. The deposit
number is what a depleted graft is short of, and leaving it out of the readout makes
the load respond to the carrier count far more steeply than it does.

## Headline numbers

- Graft load responds to carrier number as `load ~ nu^0.67`. Two independent
  perturbations give **0.70 ± 0.17** (Irf8⁻ᐟ⁻) and **0.67 ± 0.18** (BLZ945),
  weighted mean **0.687 ± 0.123** — sublinear at 2.6σ.
- **Two out-of-sample load predictions.** Irf8⁻ᐟ⁻ 0.26 against 0.244 ± 0.061 measured,
  and aged 0.26 against 0.177 ± 0.134 — neither arm fits an amyloid parameter.
  (Charging the whole of the ED Fig. 6h per-cell cargo deficit to uptake moves the
  Irf8⁻ᐟ⁻ figure to 0.15; the two readings bracket the measurement.)
- **The host load barely moves when the carriers go.** At 5% of normal carrier density
  the host load is 1.05× control while the graft falls to 0.054×. Four independent
  experiments report the host unchanged: Grathwohl et al. (2009, 95–97% ablation),
  Spangenberg et al. (2016, ~80% elimination), Olmos-Alonso et al. (2016, CSF1R
  blockade) and ED Fig. 6b of the source paper itself (41.2 vs 38.6 plaques mm⁻²).
- With carriers inert the graft load is **1.6 × 10⁻¹¹%**: the largest deposit anywhere
  inside it is 0.26 of a single seed mass, against the ten seed masses a deposit needs
  to be resolvable. The bottleneck now follows from deposit size, not from a threshold.
- Graft microglia **749** mm⁻² at 2 weeks and **518** at 4, against measured 751 ± 84
  and 517 ± 32 — the overshoot and decline are reproduced.
- Carriers import **0.24%** of the Aβ crossing into the graft at 4 weeks, **0.95%** at 16.
- Screening lengths: 0.797 mm clean, 0.406 mm graft, 0.233 mm peri-graft host.
  Export length `ell_x` = 0.186 mm; `Lambda` = 0.46.
- Peak chemotactic drift 47.6 µm day⁻¹ ≈ 2 µm h⁻¹. Nimmerjahn et al. (2005) found 5%
  of microglial somata shifting 1–2 µm h⁻¹ and no somal migration toward a lesion over
  5.5 h, and Davalos et al. (2005) found processes converging without cell-body
  movement over 1–3 h, so the fitted drift sits at the top of the measured range and is
  consistent with both null results on their timescales.
- Host load 1.86% of area at 4 weeks; graft-to-host ratio 0.31. Grathwohl et al. (2009)
  report roughly 2% Aβ-covered area in APPPS1 neocortex at the matching age.
- Lesion window amyloid rises 3.4-fold at 6 h and 4.8-fold at 24 h, against measured
  2.94 ± 1.44 and 9.40 ± 8.59.

## Known limitations of this implementation

- **The genotype ratio at the lesion is badly under-predicted**: 3.1 against a measured
  15.0 (Fig. 3c). The 24 h increments themselves are both inside the published scatter,
  which is wide, but the ratio is the tighter comparison and the model fails it.
- **Cx3cr1⁻ᐟ⁻ is not reproduced.** Those recipients have the same graft microglial
  density as controls (506 vs 517 mm⁻²) and a quarter of the amyloid (0.273, ED Fig.
  5b,c). Removing the chemotactic term entirely gives 0.58, so the model has no way to
  lose that much amyloid at an unchanged cell count.
- The host load is not measured in the source experiment, and cross-laboratory
  percentage-area comparisons are not reliable: Bader et al. (2023) give 0.08% at
  3 months for 5xFAD Aβ1-x, Grathwohl et al. (2009) about 2% for APPPS1 at 3–4 months,
  and the three cohorts within d'Errico et al. themselves span a factor of four.
- `c_sat`, the deposit density at which releases stop founding new deposits, is
  declared equal to the host density rather than measured.
- The chemoattractant amplitude is fitted separately for the injection track and for
  the laser lesion, which are different insults.
