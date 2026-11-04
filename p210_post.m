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
%     d.carrFrac   fraction of graft amyloid traceable to carrier release
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
rel = zeros(nt,1);
for k = 1:nt
    b = sol.b(k,:)';
    rel(k) = (1 - p.rho) * p.k_r * sum(G.Vc(gi) .* b(gi));
end
d.seedCum = cumtrapz(sol.t, rel);
d.carrFrac = d.seedCum ./ max(d.graftTot, eps);

% ---- mass imported across the graft boundary ----------------------------
% Inward fluxes at the face nearest r = Rg: cargo carried by cells, and
% free monomer arriving by diffusion.  Integrated in time these give the
% total mass each route has delivered into the graft.
ib = find(G.rf >= p.Rg, 1, 'first');
ib = min(max(ib,2), G.N);
qb = zeros(nt,1);  qm = zeros(nt,1);
for k = 1:nt
    m = sol.m(k,:)';  b = sol.b(k,:)';  nu = sol.nu(k,:)';  s = sol.s(k,:)';
    F = sol.F(k,:)';  ck = max(sol.c(k,:)', p.c_min);
    S = F ./ (1 + (F ./ (ck * p.m_c)).^(1/3));
    v = p.chi_s*(s(ib)-s(ib-1))/G.dr + p.chi_F*(S(ib)-S(ib-1))/G.dr;
    if v > 0, iu = ib-1; else, iu = ib; end
    fillu = max(1 - nu(iu)/p.nu_max, 0);
    Jb = -p.D_n*(b(ib)-b(ib-1))/G.dr + v*b(iu)*fillu;   % outward positive
    Jm = -p.D_m*(m(ib)-m(ib-1))/G.dr;
    qb(k) = -G.Af(ib) * Jb;                              % inward positive
    qm(k) = -G.Af(ib) * Jm;
end
d.carrImport = cumtrapz(sol.t, qb);      % cargo carried in, model units
d.diffImport = cumtrapz(sol.t, qm);      % free monomer diffused in
d.carrShare  = d.carrImport ./ max(d.carrImport + max(d.diffImport,0), eps);
d.carrPerDeposit = d.carrImport ./ max(d.graftTot, eps);

% ---- onset time against depth -------------------------------------------
d.depth  = p.Rg - G.rc(gi);
idx = find(gi);
d.onsetT = nan(numel(idx),1);
for j = 1:numel(idx)
    i = idx(j);
    % onset is the first time the deposit at this radius is resolvable and
    % covers more than p.aThr of the area
    k = find(d.areaPct(:,i) > p.aThr, 1, 'first');
    if ~isempty(k) && k > 1
        % linear interpolation in time onto the threshold crossing
        F1 = d.areaPct(k-1,i);  F2 = d.areaPct(k,i);
        d.onsetT(j) = sol.t(k-1) + (p.aThr - F1)/(F2 - F1)*(sol.t(k)-sol.t(k-1));
    elseif ~isempty(k)
        d.onsetT(j) = sol.t(k);
    end
end

% ---- transport length scales --------------------------------------------
Fend = sol.F(end,:)';  nuend = sol.nu(end,:)';
cend = max(sol.c(end,:)', p.c_min);
Send = Fend ./ (1 + (Fend ./ (cend * p.m_c)).^(1/3));
d.ell_m = sqrt(p.D_m ./ (p.lam_m + p.k_a*Send + p.k_u*nuend));
tau = 1/(p.k_d + p.k_r);
% NOTE.  These two are read off whatever output times the caller asked for,
% and the peak drift lives in the first few days, so a coarse tspan
% underestimates both badly.  Use p210_scales for any reported value.
d.Vmig = p.chi_s * max(max(abs(diff(sol.s,1,2)))) / G.dr;
d.ell_x = sqrt(2*p.D_n*tau + (d.Vmig*tau)^2);
if Vg > 0
    d.Lambda = d.ell_x / mean(d.ell_m(gi));
else
    d.Lambda = d.ell_x / mean(d.ell_m);
end
end
