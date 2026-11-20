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
meanfield = isfield(p,'readout') && strcmpi(p.readout,'meanfield');

d.areaPct = zeros(nt,G.N);
for k = 1:nt
    F = sol.F(k,:)';
    c = max(sol.c(k,:)', p.c_min);
    if meanfield
        a = F / p.perPercent;
        a(F <= p.F_thr) = 0;
        d.areaPct(k,:) = a';
        if Vg > 0
            d.loadPct(k)    = sum(G.Vc(gi) .* a(gi)) / Vg;
            d.graftTot(k)   = sum(G.Vc(gi) .* F(gi));
            d.nuGraftVol(k) = sum(G.Vc(gi) .* sol.nu(k,gi)') / Vg * p.n0_3D;
        end
        d.hostPct(k) = sum(G.Vc(hi) .* a(hi)) / Vh;
        continue
    end
    mbar = max(F ./ c, realmin);                     % mass per deposit
    % Visibility rises smoothly with deposit size rather than switching on
    % at a threshold: half the deposits of mass p.detRatio*p.m_seed are
    % called plaques.  A hard cut would put a variant on one side or the
    % other of a cliff and decide the answer by itself.
    vis = 1 ./ (1 + (p.detRatio * p.m_seed ./ mbar).^3);
    a = p.kap * c.^(1/3) .* F.^(2/3) .* vis;
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
if meanfield
    ad = Fd / p.perPercent;
    ad(Fd <= p.F_thr) = 0;
else
    ad = p.kap * cd.^(1/3) .* Fd.^(2/3) ./ (1 + (p.detRatio * p.m_seed ./ max(Fd ./ cd, realmin)).^3);
end
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
