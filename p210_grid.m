function G = p210_grid(p)
%P210_GRID  Finite-volume grid in planar, cylindrical or spherical geometry.
%
%   Cells 1..N span 0 <= r <= Rout.  Faces 1..N+1; face i separates cell i-1
%   from cell i, so face 1 is the symmetry axis (zero area in curved
%   geometry) and face N+1 is the outer no-flux boundary.
%
%   Fields: dr, rf (faces), rc (centres), Af (face areas per unit solid
%   angle), Vc (cell volumes), inGraft, inHost, Qshape, JP.

N  = p.N;
G.N  = N;
G.dr = p.Rout / N;
G.rf = (0:N)' * G.dr;
G.rc = 0.5 * (G.rf(1:N) + G.rf(2:N+1));
g    = p.geom;
G.Af = G.rf .^ g;
G.Vc = (G.rf(2:N+1).^(g+1) - G.rf(1:N).^(g+1)) / (g+1);

G.inGraft = G.rc <  p.Rg;
G.inHost  = G.rc >= p.Rg;
G.Vgraft  = sum(G.Vc(G.inGraft));
G.Vhost   = sum(G.Vc(G.inHost));

% ---- microglial carrying capacity, by compartment -----------------------
%   The graft supports fewer microglia than mature cortex; p.fnu_graft is
%   the ratio, taken from the measured graft density once the acute phase
%   has passed.  fnu_graft = 1 makes the capacity uniform.
fnu = 1; if isfield(p,'fnu_graft'), fnu = p.fnu_graft; end
G.nu0f = p.nu0 * (double(G.inHost) + fnu * double(G.inGraft));
if p.Rg <= 0, G.nu0f = p.nu0 * ones(N,1); end

% ---- shape of the injury source ----------------------------------------
if ~isfield(p,'srcmode'), p.srcmode = 'graft'; end
switch lower(p.srcmode)
    case 'graft'        % injection track and transplant trauma
        G.Qshape = double(G.inGraft);
    case 'point'        % focal laser lesion at r = 0
        w = 0.02; if isfield(p,'w_s'), w = p.w_s; end
        G.Qshape = exp(-(G.rc / w).^2);
    case 'none'
        G.Qshape = zeros(N,1);
    otherwise
        error('p210_grid:srcmode','unknown source mode ''%s''', p.srcmode);
end

% ---- Jacobian sparsity: 6 fields on a 3-point stencil -------------------
T    = spdiags(ones(N,3), -1:1, N, N);
G.JP = kron(ones(6,6), T);
end
