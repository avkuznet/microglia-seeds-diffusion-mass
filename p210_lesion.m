function L = p210_lesion(variant, Q0les, perPercent, chi_s, m_seed)
%P210_LESION  Focal laser-lesion geometry used to calibrate carrier motion.
%
%   L = P210_LESION('baseline')  or  P210_LESION('irf8ko')
%
%   Homogeneous 7-10 month 5xFAD cortex (no graft) with a point injury at
%   r = 0, imaged for 24 h.  Observables are taken inside the equivalent
%   circle of the published 150 um x 150 um window (radius 84.6 um):
%     L.tHr        times, hours
%     L.nuWin      microglia in the window, cells/mm^2 on a 25 um section
%     L.AbWin      total amyloid in the window (cargo plus deposit), model units
%     L.AbWinRel   amyloid in the window relative to its value at t = 0
%     L.respFrac   fraction of cells displaced toward the lesion, by distance bin
%     L.binEdges   distance bins used, um

if nargin < 1 || isempty(variant), variant = 'baseline'; end

p = p210_params(variant);
if nargin > 1 && ~isempty(Q0les), p.Q0les = Q0les; end
if nargin > 2 && ~isempty(perPercent), p.perPercent = perPercent; end
% The carrier field constants must be threaded in explicitly.  Reading them
% from p210_fitted.mat instead makes the lesion fit depend on whatever was
% in that file when the calibration started, so a calibration from scratch
% and one from a warm file give different answers.
if nargin > 3 && ~isempty(chi_s),  p.chi_s  = chi_s;  end
if nargin > 4 && ~isempty(m_seed), p.m_seed = m_seed; end
p.geom    = 1;            % cylindrical: a lesion imaged in a cortical plane
p.Rg      = 0;            % no graft
p.Rout    = 0.30;         % mm
p.N       = 300;
p.srcmode = 'point';
p.w_s     = 0.020;        % mm, laser spot (about 80 um^2)
p.tau_s   = 0.25;         % day, acute injury signal
p.lam_s   = 2.0;          % 1/day
p.Q0      = p.Q0les;      % 1/day, FITTED to the microglial accumulation

% Deposit levels are declared as areas and converted with perPercent.  The
% cortex around the lesion carries the established plaque load of a 7-10
% month animal; the counting window itself starts far below that average,
% at the level implied by the published amyloid area at the lesion site and
% its relative increment over 24 h (Fig. 3c with Fig. 3d).  A lesion is
% aimed at a site, not at the tissue mean, and starting the window at the
% mean is what previously made the delivered cargo invisible against it.
D = p210_data();
Rw = 0.0846;                                % equivalent radius of the window
p.F_h_pct    = 3.5;                         % surrounding cortex, % of area
p.F_hole_pct = D.lesWin0Pct;                % the window at t = 0, % of area
p            = p210_derive(p);   % F_h, F_hole, c_h, c_hole, kap, m_c
p.F_hole_r   = Rw;

tspan = linspace(0, 1, 25);                % 0 to 24 h
sol   = p210_solve(p, tspan);

G  = sol.G;
win = G.rc < Rw;
Vw  = sum(G.Vc(win));

L.tHr    = sol.t * 24;
L.nuWin  = zeros(numel(sol.t),1);
L.AbWin  = zeros(numel(sol.t),1);
for k = 1:numel(sol.t)
    L.nuWin(k) = sum(G.Vc(win) .* sol.nu(k,win)') / Vw * p.n0_3D * p.sect;
    L.AbWin(k) = sum(G.Vc(win) .* (sol.b(k,win)' + sol.F(k,win)')) / Vw;
end
L.AbWinRel = L.AbWin / L.AbWin(1);

% ---- displacement toward the lesion, binned as in the source paper ------
L.binEdges = [0 30 60 90 120];              % um
nb = numel(L.binEdges) - 1;
L.respFrac = zeros(nb,1);
n0 = sol.nu(1,:)';  nE = sol.nu(end,:)';
for j = 1:nb
    sel = G.rc*1e3 >= L.binEdges(j) & G.rc*1e3 < L.binEdges(j+1);
    if any(sel)
        % fraction of the local population that has moved, as the relative
        % change in density within the bin
        L.respFrac(j) = max(0, min(1, abs(sum(G.Vc(sel).*(nE(sel)-n0(sel)))) ...
                                       / max(sum(G.Vc(sel).*n0(sel)), eps)));
    end
end
L.p   = p;
L.sol = sol;
end
