function S = p210_scales(p, nt)
%P210_SCALES  Transport length scales, on a time grid fine enough to resolve them.
%
%   S = P210_SCALES(p)       uses 120 output times over the first 4 weeks
%   S = P210_SCALES(p, nt)   uses nt of them
%
%   Returns S.Vmig (mm/day), S.ell_x, S.ell_m_clean, S.ell_m_graft,
%   S.ell_m_host, S.Lambda, S.tau.
%
%   WHY THIS IS NOT IN p210_post.  The peak chemotactic drift occurs in the
%   first days, while the injection-injury signal is still steep, and decays
%   with it.  Reading it off the handful of output times the rest of the
%   analysis uses (0, 1, 2, 4, 16 weeks) misses the peak badly: five times
%   give 48 um/day, twenty give 91, and the value is converged by about
%   twenty.  The export length and Lambda inherit the error.  The scales are
%   therefore computed here, from a solve of their own, rather than from
%   whatever times a caller happened to ask for.

if nargin < 2 || isempty(nt), nt = 120; end

TS  = linspace(0, 28, nt);
sol = p210_solve(p, TS);
G   = sol.G;

% peak chemotactic drift over space and time
S.Vmig = p.chi_s * max(max(abs(diff(sol.s, 1, 2)))) / G.dr;

S.tau   = 1 / (p.k_d + p.k_r);
S.ell_x = sqrt(2*p.D_n*S.tau + (S.Vmig*S.tau)^2);
S.drift = S.Vmig * S.tau;              % directed contribution
S.wander = sqrt(2*p.D_n*S.tau);        % random contribution

% screening length at the final time, by compartment
F = sol.F(end,:)';  nu = sol.nu(end,:)';  c = max(sol.c(end,:)', p.c_min);
Sf = F ./ (1 + (F ./ (c * p.m_c)).^(1/3));
ellm = sqrt(p.D_m ./ (p.lam_m + p.k_a*Sf + p.k_u*nu));

S.ell_m        = ellm;
S.ell_m_clean  = sqrt(p.D_m / p.lam_m);
if G.Vgraft > 0
    S.ell_m_graft = mean(ellm(G.inGraft));
    S.Lambda      = S.ell_x / S.ell_m_graft;
else
    S.ell_m_graft = NaN;
    S.Lambda      = S.ell_x / mean(ellm);
end
shell = G.inHost & G.rc < 2*p.Rg;
S.ell_m_host = mean(ellm(shell));
end
