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
m0  = zeros(N,1);  m0(G.inHost)  = 0.4;
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
b0(G.inHost) = arrayfun(@(F) carrier_steady(p, 0.4, F, p.c_h), F0(G.inHost)) * p.nu0;
s0  = zeros(N,1);
y0  = [m0; F0; b0; nu0; s0; c0];

odeopt = odeset('RelTol',opts.RelTol, 'AbsTol',opts.AbsTol, ...
                'JPattern',G.JP, 'MaxStep', max(diff(tspan))*2 + 1);

[t, Y] = ode15s(@(t,y) p210_rhs(t, y, p, G), tspan, y0, odeopt);

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
function b = carrier_steady(p, m, F, c)
%Steady per-carrier load in tissue with free level m and deposit F.
S = F ./ (1 + (F ./ (max(c, p.c_min) * p.m_c)).^(1/3));
up = p.k_u * m + p.phi * p.k_p * S;
b  = min(up / (p.k_d + p.k_r), p.b_max);
end
