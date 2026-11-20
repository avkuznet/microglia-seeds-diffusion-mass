function sol = p210_solve(p, tspan, opts)
%P210_SOLVE  Method-of-lines solution of the carrier-diffusion model.
%
%   sol = P210_SOLVE(p, tspan)        integrates to the times in tspan (days)
%   sol = P210_SOLVE(p, tspan, opts)  with opts.RelTol, opts.AbsTol, opts.qss
%
%   State: y = [m; F; b; nu; s], each of length p.N.
%   Transport of nu and b uses the same flux-limited upwind operator and the
%   same velocity field, so that the per-carrier load b/nu cannot develop
%   spurious extrema at the migration front.
%
%   Returns sol.t, sol.m, sol.F, sol.b, sol.nu, sol.s (each nt-by-N), the grid
%   G and the parameter struct p.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'RelTol'), opts.RelTol = 1e-6; end
if ~isfield(opts,'AbsTol'), opts.AbsTol = 1e-9; end
if ~isfield(opts,'qss'),    opts.qss    = false; end

G = p210_grid(p);
N = p.N;

% ---- initial condition --------------------------------------------------
m0  = zeros(N,1);  m0(G.inHost)  = p.m_h;
F0  = zeros(N,1);  F0(G.inHost)  = p.F_h;
nu0 = zeros(N,1);  nu0(G.inHost) = p.nu0;
% A laser lesion is aimed at a site, and the window in which amyloid is
% counted starts well below the tissue average rather than at it.  Where
% p.F_hole_r is set, the deposit inside that radius starts at p.F_hole.
if isfield(p,'F_hole_r') && p.F_hole_r > 0
    F0(G.rc < p.F_hole_r) = p.F_hole;
end
% The host already carries deposits, at the declared number density p.c_h;
% the graft carries none, which is the seeding bottleneck in its sharpest
% form - not a small deposit but no deposit to grow on.
c0  = zeros(N,1);  c0(G.inHost)  = p.c_h;
if isfield(p,'F_hole_r') && p.F_hole_r > 0
    c0(G.rc < p.F_hole_r) = p.c_hole;
end
b0  = zeros(N,1);
b0(G.inHost) = arrayfun(@(F) carrier_steady(p, p.m_h, F, p.c_h), F0(G.inHost)) * p.nu0;
s0  = zeros(N,1);
y0  = [m0; F0; b0; nu0; s0; c0];

odeopt = odeset('RelTol',opts.RelTol, 'AbsTol',opts.AbsTol, ...
                'JPattern',G.JP, 'MaxStep', max(diff(tspan))*2 + 1);

%   The solver returns its own steps, and the state at the requested times
%   is read off them.  The same steps carry the time integrals of the
%   boundary fluxes and the onset times (see densediag below), which a
%   handful of output times cannot resolve.
S = ode15s(@(t,y) p210_rhs(t, y, p, G), tspan([1 end]), y0, odeopt);
t = tspan(:);
Y = deval(S, t')';

sol.dense = densediag(S, p, G, t);
sol.t  = t;
sol.m  = Y(:,      1:N  );
sol.F  = Y(:,   N+1:2*N );
sol.b  = Y(:, 2*N+1:3*N );
sol.nu = Y(:, 3*N+1:4*N );
sol.s  = Y(:, 4*N+1:5*N );
sol.c  = Y(:, 5*N+1:6*N );
sol.G  = G;
sol.p  = p;
end

% -------------------------------------------------------------------------
function dn = densediag(S, p, G, t)
%Quantities that need the whole time course, not the state at a few times.
%
%   The fluxes across the graft boundary are integrated over the solver's
%   own steps, with three Gauss-Legendre nodes in each.  The diffusive flux
%   is very large for a short time after t = 0, when free amyloid first
%   meets the empty graft, and the solver takes small steps there; a
%   trapezoid rule on output times a week apart overstates the diffusive
%   import many times over.  The deposit fields in the graft are also kept
%   on a fine uniform time grid, for the onset times.
%
%   dn.carrImport, dn.diffImport   cumulative inward transport of cargo and
%                                  of free amyloid across r = Rg, at times t
%   dn.seedCum                     cumulative insoluble carrier release
%                                  inside the graft, at times t
%   dn.tt, dn.Fg, dn.cg            fine time grid; F and c in the graft cells
%   dn.gradsMax                    largest gradient of the injury signal over
%                                  space and over the solver's steps; the drift
%                                  peaks in the first days, between output times
N  = G.N;
gi = find(G.inGraft);  gi = gi(:)';
%   The boundary is the face between the last graft cell and the first host
%   cell.  (The first face at or beyond r = Rg is that face only when Rg/dr
%   is a whole number or rounds up; otherwise it lies one cell out in the
%   host.)
ib = numel(gi) + 1;
ib = min(max(ib,2), N);

tb = unique([S.x(:); t(:)]);
tb = tb(tb >= t(1) & tb <= t(end));
h  = diff(tb);
g  = [-sqrt(3/5), 0, sqrt(3/5)];  w = [5; 8; 5] / 9;
tq = (tb(1:end-1) + (h/2) * (1 + g))';        % 3 nodes per interval, in columns
tq = tq(:)';

cl = [ib-1, ib];
Yq = deval(S, tq, [cl, N+cl, 2*N+cl, 3*N+cl, 4*N+cl, 5*N+cl, 2*N+gi]);
m = Yq(1:2,:);   F = Yq(3:4,:);   b = Yq(5:6,:);
nu = Yq(7:8,:);  s = Yq(9:10,:);  c = max(Yq(11:12,:), p.c_min);
Sf = F ./ (1 + (F ./ (c * p.m_c)).^(1/3));
v  = p.chi_s*(s(2,:) - s(1,:))/G.dr + p.chi_F*(Sf(2,:) - Sf(1,:))/G.dr;
up = v > 0;                                    % outward flow: donor is the inner cell
nuu = nu(2,:);  nuu(up) = nu(1,up);
bu  = b(2,:);   bu(up)  = b(1,up);
Jb = -p.D_n*(b(2,:) - b(1,:))/G.dr + v .* bu .* max(1 - nuu/p.nu_max, 0);   % outward positive
Jm = -p.D_m*(m(2,:) - m(1,:))/G.dr;
qb = -G.Af(ib) * Jb;                           % inward positive
qm = -G.Af(ib) * Jm;
rel = (1 - p.rho) * p.k_r * (G.Vc(gi)' * Yq(13:end,:));

cum = @(f) [0; cumsum((h/2) .* (reshape(f, 3, [])' * w))];
[~, loc] = ismember(t, tb);
Ib = cum(qb);  Im = cum(qm);  Ir = cum(rel);
dn.carrImport = Ib(loc);
dn.diffImport = Im(loc);
dn.seedCum    = Ir(loc);

inT = S.x >= t(1) & S.x <= t(end);
dn.gradsMax = max(max(abs(diff(S.y(4*N+1:5*N, inT), 1, 1)))) / G.dr;

nf = min(max(ceil((t(end) - t(1)) / 0.02), 1), 20000);
dn.tt = linspace(t(1), t(end), nf + 1)';
Yg = deval(S, dn.tt', [N+gi, 5*N+gi]);
dn.Fg = Yg(1:numel(gi), :)';
dn.cg = Yg(numel(gi)+1:end, :)';
end

% -------------------------------------------------------------------------
function b = carrier_steady(p, m, F, c)
%Steady per-carrier load in tissue with free level m and deposit F.
S = F ./ (1 + (F ./ (max(c, p.c_min) * p.m_c)).^(1/3));
up = p.k_u * m + p.phi * p.k_p * S;
b  = min(up / (p.k_d + p.k_r), p.b_max);
end
