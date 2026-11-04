function p = p210_params(variant)
%P210_PARAMS  Parameter sets for the microglial carrier-diffusion model.
%
%   p = P210_PARAMS            returns the baseline (young-adult recipient) set
%   p = P210_PARAMS(variant)   returns a perturbed set.  Variants:
%       'baseline'   young-adult Cx3cr1+/- 5xFAD recipient
%       'irf8ko'     Irf8-/- recipient: impaired somal chemotaxis and motility
%       'depleted'   BLZ945 CSF-1R inhibition: 80% fewer microglia
%       'aged'       8-11 month recipient: reduced phagocytic competence
%       'nocarrier'  carriers removed entirely (diffusion-only control)
%
%   UNITS.  Length mm, time days.  All amyloid-beta variables (m, F, b) are
%   expressed in a single common unit: the concentration of free Abeta in
%   clean tissue, m_ref.  The conversion to the histological readout is
%   p.perPercent units of F per 1% 6E10-positive area fraction, so that
%   observable load = F / p.perPercent.  Microglial density nu is expressed in
%   multiples of the homeostatic cortical density p.n0_3D.
%
%   Parameters flagged FITTED are set by p210_calibrate and stored in
%   p210_fitted.mat; everything else is measured, derived, or declared.
%
%   A. V. Kuznetsov, NC State University.  MIT licence.

if nargin < 1 || isempty(variant), variant = 'baseline'; end
p.variant = lower(variant);

% ---- geometry -----------------------------------------------------------
p.geom = 2;            % 0 planar, 1 cylindrical, 2 spherical
p.Rg   = 0.45;         % graft radius, mm
p.Rout = 2.50;         % outer radius of the simulated host shell, mm
p.N    = 320;          % finite-volume cells

% ---- free amyloid-beta --------------------------------------------------
p.D_m   = 5.383;       % mm^2/day  (62.3 um^2/s, Waters 2010)
p.lam_m = 8.467;       % 1/day     (9.8e-5 1/s, t_half 1.97 h, Cirrito 2003)
p.q0    = 8.467;       % host production; normalises m -> 1 in clean tissue
p.deplTarget = 0.80;   % target host monomer depletion at 4 weeks (Hong 2011)
p.rho   = 0.10;        % soluble fraction of released carrier cargo
p.k_1   = 1e-8;        % homogeneous nucleation, per day per m^2
p.k_2   = 0.50;        % secondary nucleation on accessible surface
p.S_c   = Inf;         % sub-voxel competition scale, model units of surface.
                       % Growth stops being proportional to the surface
                       % present once the deposits sharing a voxel consume
                       % monomer as fast as it is resupplied: consumption is
                       % k_a m S per unit volume against resupply lam_m m, so
                       % the two balance near S_c = lam_m/k_a ~ 30 with the
                       % measured clearance rate.  Declared, not fitted,
                       % and not swept in any published figure; the paper
                       % discusses it in the limitations section as a
                       % model-structure sensitivity, and p210_competition
                       % is where the alternative closure is exercised.
                       % S_c = Inf is the mean-field closure of the first
                       % version of this model.

% ---- deposited amyloid --------------------------------------------------
%   DEPOSIT NUMBER.  Amyloid is laid down as discrete deposits and the
%   histological readout is the area they cover, not the mass they contain.
%   The model therefore carries their number density c as a sixth field.
%   The observation operator is Eq. (16'):  area = kap * c^(1/3) * F^(2/3),
%   which is c deposits each of mass F/c, each presenting an area
%   proportional to its mass to the two-thirds.  A deposit is counted only
%   if it has grown to p.detRatio times the seed it started from.
p.F_hRef_pct = 0.20;   % reference host load that fixes kap, % of area
p.c_hRef = 857;        % reference host deposit density that fixes kap
p.c_h   = 857;         % MEASURED  deposit number density in the host, per
                       % mm^3.  ED Fig. 6b counts 38.6 +/- 5.7 cortical Abeta
                       % plaques per mm^2 of section; with a 25 um section
                       % and a 20 um plaque the Abercrombie correction gives
                       % 38.6/(0.025+0.020) = 857 per mm^3.  The same figure
                       % finds 41.2 +/- 5.9 in Irf8-/- recipients, so the
                       % host deposit density does not depend on the carrier
                       % population - only the naive graft does.
p.c_sat = 857;         % DECLARED  density at which new deposits stop being
                       % made because a release lands on one already there.
                       % Equal to p.c_h: the host has reached its own
                       % saturation density, which is why it has that one.
p.c_min = 1.0e-12;     % floor, so that S and the readout stay finite
p.m_seed = 9.375e-4;   % FITTED  deposit mass left by one carrier release
                       % event, model units per deposit; fitted to the
                       % BLZ945 load ratio.  It has an independent estimate
                       % it is NOT fitted to: a loaded microglial cell
                       % carries about one Abeta cluster (ED Fig. 4d, 1.04
                       % per cell) and leaves it behind when it releases, so
                       % one release event should start one deposit of mass
                       % p.b_max/p.n0_3D = 9.375e-4.  The fitted value lands
                       % at 0.49 of that, which is the closest thing in this
                       % model to a check on the seeding picture itself.
p.seedRatio = 1.0e3;   % DECLARED  centres made per unit mass by homogeneous
                       % nucleation, relative to carrier release.  A nucleus
                       % is a few monomers and a released cluster is many
                       % thousands, so homogeneous nucleation makes far more
                       % centres, each of them far too small to resolve.
p.detRatio = 10;       % DECLARED  growth factor over the seed mass at which
                       % a deposit becomes resolvable
p.kap   = 0.02;        % DERIVED below from c_h and F_h
p.perPercent = 500;    % FITTED  model units of F per 1% area; the single
                       % unit conversion, fitted to the 4-week graft load
p.m_c   = 2.5;         % DERIVED below, p.F_c / p.c_h
p.F_c   = 500;         % compaction scale of a deposit population at c_h
p.F_h_pct = 0.20;      % host load at transplantation, % of area.  9-10 week
                       % 5xFAD recipients are pre-depositing; deposit levels
                       % are declared as areas and converted with perPercent,
                       % so that no deposit is specified in model units alone.
p.F_h   = 100;         % host load in model units (set from F_h_pct below)
p.F_thr = 25;          % mean-field detection threshold, retained for the
                       % mean-field comparison in p210_competition
p.aThr  = 0.05;        % smallest load called detectable, % of area
p.k_c   = 0;           % extra non-phagocytic deposit clearance (phagocytosis is explicit)

% ---- carriers -----------------------------------------------------------
p.nu_max = 4.0;        % volume-filling ceiling, multiples of n0
p.nu0    = 1.0;        % homeostatic density in the host, multiples of n0
p.D_n    = 2.246e-3;   % mm^2/day random motility (0.026 um^2/s)
p.chi_s  = 5.8e-2;     % FITTED  chemotaxis to injury signal, mm^2/day per
                       % unit s, from the 2-week graft density (Fig. 1d)
p.chi_F  = 0.0;        % mm^2/day per unit S, recruitment to deposit surface
p.g_n    = 0.05;       % 1/day relaxation of density toward the local
                       % carrying capacity
p.fnu_graft = 0.35;    % FITTED  microglial carrying capacity of the graft,
                       % as a fraction of the homeostatic cortical density,
                       % from the 4-week graft density (Fig. 1d).  The graft
                       % is embryonic tissue and supports fewer microglia.
p.del_n  = 0.0;        % 1/day death (TUNEL-negative: set to zero)
p.k_u    = 0.10;       % uptake of soluble Abeta, 1/day per nu
p.k_p    = 0.040;      % phagocytosis of deposit surface, 1/day per nu per S
p.b_max  = 30;         % cargo capacity per unit nu
p.k_d    = 0.030;      % lysosomal degradation of cargo, 1/day
p.k_r    = 0.300;      % release of cargo, 1/day
p.phi    = 1.0;        % phagocytic competence of the carrier population
p.phiAged = 0.494;     % EX VIVO PROXY, not fitted and not a direct uptake
                       % measurement: the ratio of intracellular Abeta
                       % cluster burden in old against adult isolated
                       % microglia, ED Fig. 4d.  ED Fig. 4a-d is confocal
                       % microscopy and cell counting after isolation and
                       % 24 h culture (panels e-g are the flow-cytometry
                       % series), so it records cargo present at the end of
                       % culture rather than an uptake flux.  It also comes
                       % from mice of both sexes, unlike the in vivo
                       % grafting experiments, and from an older cohort
                       % (50-60 weeks against 20-30) than the 8-11 month
                       % recipients it is applied to.
p.krAged  = 0.0207;    % NOT USED by default.  The aged release rate the
                       % previous version fitted; see the 'aged' variant.

% ---- Irf8-/- recipients, all three measured -----------------------------
p.nu0IRF8    = 0.133;  % MEASURED  graft microglial density as a fraction of
                       % control: 68.9 against 519.5 per mm^2, Fig. 2g
p.cargoIRF8  = 0.663;  % MEASURED  amyloid carried per cell, as a fraction
                       % of control: 4.75 against 7.17, ED Fig. 6h

% ---- injury signal ------------------------------------------------------
p.D_s   = 0.05;        % mm^2/day
p.lam_s = 0.30;        % 1/day
p.Q0    = 0.60;        % 1/day, source strength over the graft
p.tau_s = 7.0;         % day, decay of the injection-injury source

% ---- observation --------------------------------------------------------
p.sect  = 0.025;       % histological section thickness, mm
p.n0_3D = 3.2e4;       % homeostatic microglial density, cells/mm^3

% ---- fitted values, if a calibration has been run -----------------------
p.f_chi = 0.15;        % FITTED  Irf8-/- multiplier on chemotaxis and motility
p.k_a   = 0.40;        % FITTED  surface growth rate constant, 1/day per unit S
p.Q0les = 600;         % FITTED  acute lesion signal strength, 1/day
if exist('p210_fitted.mat','file')
    L = load('p210_fitted.mat');
    fn = {'k_a','f_chi','Q0les','perPercent','F_thr', ...
          'chi_s','fnu_graft','m_seed'};
    for q = 1:numel(fn)
        if isfield(L,fn{q}), p.(fn{q}) = L.(fn{q}); end
    end
end

% ---- variants -----------------------------------------------------------
switch p.variant
    case 'baseline'
        % nothing
    case 'irf8ko'
        % The carrier field is imposed at its measured value and the amyloid
        % load is predicted, exactly as for BLZ945.  Whether Irf8-/- recipients
        % have fewer cortical microglia or cells that fail to enter the graft
        % is unsettled - Horiuchi et al. (2012) report colonisation to be
        % largely normal - and the model does not need to know: what it needs
        % is how many carriers end up in the graft, which Fig. 2g measures.
        %
        % Motility is NOT scaled.  Fig. 2d,e measure the travel of microglial
        % PROCESSES, and the model's drift is the translocation of the SOMA.
        % Davalos et al. (2005) and Nimmerjahn et al. (2005) show the two are
        % different quantities: processes converge on a lesion at ~1.5 um/min
        % while the cell bodies do not move at all over several hours.
        %
        % ED Fig. 6h finds 0.663 as much amyloid in each Irf8-/- cortical
        % microglial cell.  That is a standing content, k_u/(k_d+k_r), not an
        % uptake rate, and the model cannot tell from it which of the three
        % rates moved: a cell that releases more slowly would show more per
        % cell, not less.  Charging the whole deficit to uptake, which
        % p.useIrf8Cargo does, takes the predicted load ratio from 0.256 to
        % 0.149 against a measured 0.244 +/- 0.061, so the two readings
        % bracket the measurement.  The default charges it to nothing, which
        % makes this arm a pure carrier-number perturbation like BLZ945 -
        % and the two are measured to have the same cell ratio, 0.133 and
        % 0.132, and the same load ratio, 0.244 and 0.256.
        p.nu0 = p.nu0IRF8;                % MEASURED  Fig. 2g
        if isfield(p,'useIrf8Cargo') && p.useIrf8Cargo
            p.k_u = p.cargoIRF8 * p.k_u;
            p.phi = p.cargoIRF8 * p.phi;
        end
    case 'depleted'
        % graft microglial density measured after BLZ945, ED Fig. 7f:
        % 77.4 against 585.9 per mm^2, a 87% reduction rather than the
        % 80% quoted for the cortex as a whole.
        p.nu0 = 0.132;
    case 'aged'
        % NOTHING IS FITTED HERE.  Uptake per cell is set to the ex vivo
        % proxy of 0.494 (ED Fig. 4d, a cluster-burden ratio rather than a
        % measured uptake rate) and the host load is the declared one for an
        % 8-11 month recipient; the graft load is then a prediction,
        % conditional on that proxy.
        %
        % Earlier versions of this model also fitted a reduced cargo release
        % rate for aged carriers, because with a mass-proportional readout
        % the measured two-fold uptake deficit could not come close to the
        % observed load ratio - it needed a thirty-fold loss.  With the
        % discrete-deposit readout the load follows delivery sublinearly,
        % the measured deficit alone lands inside the experimental error,
        % and the extra fitted constant is no longer needed.  Setting
        % p.useAgedRelease applies it, which moves the predicted ratio from
        % 0.254 to 0.175 against a measured 0.177 +/- 0.134.
        p.phi     = p.phiAged;            % EX VIVO PROXY  ED Fig. 4d
        p.F_h_pct = 4.0;                  % 8-11 month recipients, % area
        if isfield(p,'useAgedRelease') && p.useAgedRelease
            p.k_r = p.krAged;
        end
    case 'cx3cr1ko'
        % Cx3cr1-/- recipients have a normal graft microglial density, 506
        % against 517 per mm^2 (ED Fig. 5b), and a quarter of the amyloid
        % (ED Fig. 5c).  The model already distinguishes the two routes by
        % which a carrier loads: U_m from the soluble pool, at k_u, and U_F
        % from the deposit surface, at phi k_p.  Cx3cr1 mediates selective
        % phagocytosis of protofibrillar Abeta (Liu et al. 2010), which is
        % the second route, so the genotype is represented by switching that
        % route off and changing nothing else.  phi = 0 is a limit, not a
        % fitted value, and the carrier density is left alone.
        p.phi = 0;
    case 'nocarrier'
        p.k_p = 0; p.k_u = 0; p.k_r = 0;  % carriers present but inert
    otherwise
        error('p210_params:variant','unknown variant ''%s''', variant);
end

% ---- deposit levels declared as areas, converted to model units ---------
p = p210_derive(p);
end
