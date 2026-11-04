function dy = p210_rhs(t, y, p, G)
%P210_RHS  Right-hand side of the carrier-diffusion system.
%
%   Amyloid-beta is conserved between the three pools m, F and b apart from
%   production q0 in the host, first-order clearance of the free pool, and
%   lysosomal degradation of cargo.  Material taken up by a carrier is
%   removed from the pool it came from; material released by a carrier is
%   returned, a fraction rho of it soluble.

N  = G.N;
dr = G.dr;
Af = G.Af;
Vc = G.Vc;

m  = max(y(      1:N  ), 0);
F  = max(y(   N+1:2*N ), 0);
b  = max(y( 2*N+1:3*N ), 0);
nu = max(y( 3*N+1:4*N ), 0);
s  = max(y( 4*N+1:5*N ), 0);
c  = max(y( 5*N+1:6*N ), 0);      % number density of distinguishable deposits

% ---- accessible fibril surface ------------------------------------------
%   A deposit of mass M exposes a surface M/(1+(M/m_c)^(1/3)): proportional
%   to its mass while it is small, to M^(2/3) once it is compact.  With c
%   deposits per unit volume each of mass F/c, the surface per unit volume
%   is c times that, which is the expression below.  The mean-field form of
%   the first version of this model is the special case c = p.c_h, with
%   p.F_c = p.c_h * p.m_c: it holds the deposit number fixed at the host
%   value and therefore cannot see what a carrier perturbation does to it.
S = F ./ (1 + (F ./ (max(c, p.c_min) * p.m_c)).^(1/3));

% ---- injury signal source ------------------------------------------------
Q = p.Q0 * G.Qshape * exp(-t / p.tau_s);

% ---- diffusion of m and s ------------------------------------------------
Lm = fvlap(m, p.D_m, Af, Vc, dr, N);
Ls = fvlap(s, p.D_s, Af, Vc, dr, N);

% ---- chemotactic velocity at interior faces ------------------------------
vf      = zeros(N+1,1);
vf(2:N) = p.chi_s * (s(2:N) - s(1:N-1)) / dr ...
        + p.chi_F * (S(2:N) - S(1:N-1)) / dr;

% ---- carrier and cargo transport -----------------------------------------
[Ldn, Ldb] = celltrans(nu, b, vf, p.D_n, p.nu_max, Af, Vc, dr, N);

% ---- deposition: elongation plus secondary nucleation on the surface -----
%   Sub-voxel competition.  The mean-field form (k_a m + k_2 m^2) S makes the
%   growth rate strictly proportional to the surface present, so a crowded
%   field of deposits grows as fast per unit surface as an isolated one.  It
%   cannot: the deposits in a voxel draw on one shared monomer supply, and
%   once the surface in a voxel is large enough the growth is limited by how
%   fast monomer reaches it rather than by how much surface there is.  The
%   factor 1/(1 + S/S_c) interpolates between the two, with S_c the surface
%   at which local consumption matches local resupply.  S_c = Inf recovers
%   the mean-field closure.
grow = (p.k_a * m + p.k_2 * m.^2) .* S ./ (1 + S / p.S_c);

% ---- carrier uptake, capped by the per-cell cargo ceiling ----------------
beta = b ./ max(nu, 1e-12);
fill = max(1 - beta / p.b_max, 0);
U_m  = p.k_u * nu .* m .* fill;                 % from the free pool
U_F  = p.phi * p.k_p * nu .* S .* fill;         % from the deposit
rel  = p.k_r * b;

% ---- balances -------------------------------------------------------------
dm  = Lm + p.q0 * G.inHost - p.lam_m * m - grow - U_m + p.rho * rel;
dF  = grow + (1 - p.rho) * rel + p.k_1 * m.^2 - U_F - p.k_c * nu .* S;
db  = Ldb + U_m + U_F - (p.k_d + p.k_r) * b;
% Microglial density relaxes toward the local carrying capacity.  The graft
% is embryonic tissue and supports fewer microglia than mature cortex: the
% measured graft density settles at 0.65 of the homeostatic cortical value
% while the acute injury transiently drives it above that.  A single
% capacity for both compartments can only rise monotonically to it.
dnu = Ldn + p.g_n * nu .* (1 - nu ./ G.nu0f) - p.del_n * nu;
ds  = Ls + Q - p.lam_s * s;

% ---- number density of deposits -----------------------------------------
%   Amyloid is not laid down as a continuum: it is laid down as discrete
%   deposits, and the histological readout is the area those deposits cover,
%   not the mass they contain.  The two differ whenever the number of
%   deposits changes, which is exactly what a carrier perturbation does.  A
%   carrier that releases its cargo leaves behind one aggregate of mass
%   p.m_seed - the Abeta cluster counted in ED Fig. 4d, of which a loaded
%   cell carries about one.  Homogeneous nucleation instead makes nuclei of
%   a few monomers, so per unit mass nucleated it makes p.seedRatio times as
%   many centres, each of them far too small to resolve.  Deposits grow but
%   are immobile, so c has no transport term.
%   A release that lands where a deposit already sits adds to it instead of
%   starting another, so the rate of new centres is damped as the local
%   density approaches saturation.  That is why the host sits at c_sat.
room = max(1 - c / p.c_sat, 0);
dc  = room .* ( (1 - p.rho) * p.k_r * b ...
              + p.seedRatio * p.k_1 * m.^2 ) / p.m_seed;

dy = [dm; dF; db; dnu; ds; dc];
end

% -------------------------------------------------------------------------
function L = fvlap(u, D, Af, Vc, dr, N)
%Conservative finite-volume Laplacian; zero flux at both boundaries.
flux      = zeros(N+1,1);
flux(2:N) = -D * (u(2:N) - u(1:N-1)) / dr;
L = -(Af(2:N+1).*flux(2:N+1) - Af(1:N).*flux(1:N)) ./ Vc;
end

% -------------------------------------------------------------------------
function [Ln, Lb] = celltrans(nu, b, vf, Dn, nu_max, Af, Vc, dr, N)
%Random motility plus volume-filling chemotaxis, upwinded.  The cargo b is
%advected by the identical operator so that b/nu stays bounded.
Jn = zeros(N+1,1);
Jb = zeros(N+1,1);

Jn(2:N) = -Dn * (nu(2:N) - nu(1:N-1)) / dr;
Jb(2:N) = -Dn * ( b(2:N) -  b(1:N-1)) / dr;

v   = vf(2:N);
up  = v > 0;                              % outward flow: donor is the inner cell
iL  = (1:N-1)';  iR = (2:N)';
idx = iR;  idx(up) = iL(up);
fill = max(1 - nu(idx) / nu_max, 0);
Jn(2:N) = Jn(2:N) + v .* nu(idx) .* fill;
Jb(2:N) = Jb(2:N) + v .*  b(idx) .* fill;

Ln = -(Af(2:N+1).*Jn(2:N+1) - Af(1:N).*Jn(1:N)) ./ Vc;
Lb = -(Af(2:N+1).*Jb(2:N+1) - Af(1:N).*Jb(1:N)) ./ Vc;
end
