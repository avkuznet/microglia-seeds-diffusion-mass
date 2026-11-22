function d = p210_post(sol)
%P210_POST  Diagnostics for a carrier-diffusion solution.
%
%   d = P210_POST(sol) returns
%     d.t          times (days)
%     d.loadPct    graft amyloid load, % area, observation operator Eq. (16)
%     d.hostPct    host amyloid load, % area
%     d.nuGraft    mean graft microglial density, cells/mm^2 on 25 um sections
%     d.nuGraftVol mean graft microglial density, cells/mm^3
%     d.seedCum    cumulative carrier-released amyloid in the graft
%     d.graftTot   total graft amyloid (model units, volume integral)
%     d.carrFrac   cumulative insoluble carrier release in the graft divided by
%                  the deposited graft mass (a diagnostic ratio, not a
%                  fraction of deposits or of their molecules)
%     d.onsetT     onset time (days) at each radius; NaN where never detected
%     d.depth      depth into the graft, Rg - r, for the onset curve
%     d.ell_m      screening length profile at the final time, mm
%     d.ell_x      export length, mm
%     d.Lambda     ell_x / ell_m evaluated in the graft
%     d.Vmig       representative carrier migration speed, mm/day

p = sol.p;  G = sol.G;
nt = numel(sol.t);
d.t = sol.t;

gi = G.inGraft;  hi = G.inHost;
Vg = G.Vgraft;   Vh = G.Vhost;

% ---- observation operator, Eq. (16) -------------------------------------
d.loadPct = zeros(nt,1);  d.hostPct = zeros(nt,1);
d.graftTot = zeros(nt,1); d.nuGraftVol = zeros(nt,1);
%   Eq. (16'):  c deposits per unit volume, each of mass F/c, each covering
%   an area proportional to its mass to the two-thirds, so the area fraction
%   is kap * c^(1/3) * F^(2/3).  A deposit is resolved only once it has
%   grown to p.detRatio times the seed it started from; mass spread over
%   centres too small to see contributes nothing, which is what keeps the
%   trace of homogeneously nucleated amyloid out of the readout.
%
%   p.readout = 'meanfield' restores the operator of the first version of
%   this model - area proportional to the deposit mass, gated by a hard
%   threshold on it, with the deposit number never entering - so that the
%   two can be compared on the same solution.
%
%   Three further switches, all off by default, are used by p210_readouts to
%   show how far the results depend on this rule (see obsArea below):
%     p.readout = 'section'   the section formula with thickness p.sectT, of
%                             which Eq. (14) is the thick-section limit and
%                             area proportional to mass the thin-section one
%     p.visExp                exponent of the visibility factor (3 in the paper)
%     p.cargoRead             1 or 2: intracellular cargo is read as well
d.areaPct = zeros(nt,G.N);
for k = 1:nt
    F = sol.F(k,:)';
    c = max(sol.c(k,:)', p.c_min);
    a = obsArea(p, F, c, sol.b(k,:)', sol.nu(k,:)');
    d.areaPct(k,:) = a';
    if Vg > 0
        d.loadPct(k)    = sum(G.Vc(gi) .* a(gi)) / Vg;
        d.graftTot(k)   = sum(G.Vc(gi) .* F(gi));
        d.nuGraftVol(k) = sum(G.Vc(gi) .* sol.nu(k,gi)') / Vg * p.n0_3D;
    end
    d.hostPct(k) = sum(G.Vc(hi) .* a(hi)) / Vh;
end
d.nuGraft = d.nuGraftVol * p.sect;        % cells per mm^2 on a 25 um section

% ---- carrier-released amyloid delivered into the graft ------------------
%   The time integral is taken by p210_solve over the solver's own steps
%   (sol.dense), so it does not depend on the output times asked for.
dn = sol.dense;
d.seedCum = dn.seedCum;
d.carrFrac = d.seedCum ./ max(d.graftTot, eps);

% ---- mass imported across the graft boundary ----------------------------
% Inward fluxes at the face nearest r = Rg: cargo carried by cells, and
% free monomer arriving by diffusion.  Integrated in time these give the
% total mass each route has delivered into the graft.  The diffusive flux is
% very large just after t = 0, so the integrals come from p210_solve too.
d.carrImport = dn.carrImport;            % cargo carried in, model units
d.diffImport = dn.diffImport;            % free monomer diffused in
d.carrShare  = d.carrImport ./ max(d.carrImport + max(d.diffImport,0), eps);
d.carrPerDeposit = d.carrImport ./ max(d.graftTot, eps);

% ---- onset time against depth -------------------------------------------
%   Read off the fine time grid kept by p210_solve, with the same
%   observation operator as above.
d.depth  = p.Rg - G.rc(gi);
idx = find(gi);
d.onsetT = nan(numel(idx),1);
Fd = max(dn.Fg, 0);  cd = max(dn.cg, p.c_min);
%   The fine time grid keeps F and c only, so the onset is that of the
%   deposits; cargo, if it is read at all (p.cargoRead), is left out here.
ad = obsArea(p, Fd, cd, [], []);
for j = 1:numel(idx)
    % onset is the first time the deposit at this radius is resolvable and
    % covers more than p.aThr of the area
    k = find(ad(:,j) > p.aThr, 1, 'first');
    if ~isempty(k) && k > 1
        % linear interpolation in time onto the threshold crossing
        F1 = ad(k-1,j);  F2 = ad(k,j);
        d.onsetT(j) = dn.tt(k-1) + (p.aThr - F1)/(F2 - F1)*(dn.tt(k)-dn.tt(k-1));
    elseif ~isempty(k)
        d.onsetT(j) = dn.tt(k);
    end
end

% ---- transport length scales --------------------------------------------
Fend = sol.F(end,:)';  nuend = sol.nu(end,:)';
cend = max(sol.c(end,:)', p.c_min);
Send = Fend ./ (1 + (Fend ./ (cend * p.m_c)).^(1/3));
d.ell_m = sqrt(p.D_m ./ (p.lam_m + p.k_a*Send + p.k_u*nuend));
tau = 1/(p.k_d + p.k_r);
% The peak drift lives in the first few days, between the usual output
% times, so it is taken from the solver's own steps (sol.dense), over the
% whole run.  p210_scales reports the same quantity over the first 4 weeks.
d.Vmig = p.chi_s * dn.gradsMax;
d.ell_x = sqrt(2*p.D_n*tau + (d.Vmig*tau)^2);
if Vg > 0
    d.Lambda = d.ell_x / mean(d.ell_m(gi));
else
    d.Lambda = d.ell_x / mean(d.ell_m);
end
end

% -------------------------------------------------------------------------
function a = obsArea(p, F, c, b, nu)
%OBSAREA  Stained area, % of section, from the deposit mass F and number c.
%
%   'discrete'   Eq. (14) of the paper: kap * c^(1/3) * F^(2/3) * vis.
%   'section'    c deposits per unit volume, each a sphere of diameter d,
%                projected through a section of thickness T:
%                    area fraction = c * (pi/4) d^2 * (T + 2d/3) * vis
%                The first term counts the spheres whose centre lies in the
%                section and the second the caps of those cut by its faces.
%                For d << T it is Eq. (14); for T = 0 it is the volume
%                fraction, proportional to mass, which is what a true planar
%                section measures.  d follows from the deposit mass as
%                d_ref * (mass / mu_ref)^(1/3), with d_ref fixed by the
%                reference host in p210_derive, so nothing new is assumed.
%   'meanfield'  area proportional to mass, hard threshold, no deposit
%                number: the operator of the first version of this model.
%
%   b and nu, when given and p.cargoRead > 0, add the intracellular cargo as
%   one object per cell: through the same visibility rule (1) or all of it (2).
ro = 'discrete';  if isfield(p,'readout'), ro = lower(p.readout); end
if strcmp(ro, 'meanfield')
    a = F / p.perPercent;
    a(F <= p.F_thr) = 0;
    return
end
ne = 3;  if isfield(p,'visExp'), ne = p.visExp; end
mbar = max(F ./ c, realmin);                     % mass per deposit
% Visibility rises smoothly with deposit size rather than switching on at a
% threshold: half the deposits of mass p.detRatio*p.m_seed are called
% plaques.  A hard cut would put a variant on one side or the other of a
% cliff and decide the answer by itself.
vis = 1 ./ (1 + (p.detRatio * p.m_seed ./ mbar).^ne);
switch ro
    case 'discrete'
        a = p.kap * c.^(1/3) .* F.^(2/3) .* vis;
    case 'section'
        dd = p.d_ref * (mbar / p.mu_ref).^(1/3);           % diameter, mm
        a  = 100 * c .* (pi/4) .* dd.^2 .* (p.sectT + 2*dd/3) .* vis;
        a(F <= 0) = 0;
    otherwise
        error('p210_post:readout', 'unknown readout ''%s''', ro);
end
if isfield(p,'cargoRead') && p.cargoRead > 0 && ~isempty(b)
    n  = max(nu * p.n0_3D, realmin);  b = max(b, 0);       % cells per mm^3
    mb = max(b ./ n, realmin);                             % cargo per cell
    if p.cargoRead == 1
        vb = 1 ./ (1 + (p.detRatio * p.m_seed ./ mb).^ne);
    else
        vb = 1;
    end
    a = a + p.kap * n.^(1/3) .* b.^(2/3) .* vb;
end
end
