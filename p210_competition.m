function C = p210_competition()
%P210_COMPETITION  The two observation operators, side by side.
%
%   C = P210_COMPETITION
%
%   SUPERSEDED by p210_readouts, which is the comparison the paper reports
%   (Table S3).  This routine is kept for the record and is NOT an equally
%   calibrated comparison; see below.
%
%   The histological readout is the area that amyloid deposits cover.  The
%   first version of this model took that area to be proportional to the
%   deposit MASS in a voxel, gated by a threshold on it, which is a mean
%   field over a deposit population whose number never enters.
%
%   The model now carries the deposit number density c and reads
%
%       area  =  kap * c^(1/3) * F^(2/3) * vis(F/c)
%
%   which is c deposits of mass F/c, each showing an area proportional to
%   its mass to the two-thirds, each visible once it has grown past a few
%   times the seed it started from.  The mass-proportional operator is the
%   special case in which the deposit number is pinned at the host value.
%
%   This routine refits TWO constants for each operator - k_a to the host
%   monomer depletion and perPercent to the 4-week graft load - and reports
%   what each then gives.  It does NOT refit m_seed: that constant was
%   fitted to the BLZ945 load ratio under the discrete operator and is used
%   unchanged under the mean-field one.  The BLZ945 ratio is therefore
%   fitted in the left column and predicted in the right, and the two
%   operators also differ in how a deposit is detected (a smooth visibility
%   factor against a hard threshold).  The columns are not like for like,
%   and no conclusion about which operator performs better should be drawn
%   from them.  p210_readouts refits the same three targets under each
%   rule, with the same visibility factor.

D  = p210_data();
TS = [0 7 14 21 28];
T16 = [0 7 14 28 112];
NC = 400;                % the reference grid of p210_params
TGT = D.irf8Load.val(1);

modes = {'discrete','meanfield'};
C = struct();

fprintf('\n=== observation operator: deposit number in, or left out ===\n\n');
fprintf('%-26s %10s %10s\n', '', 'discrete', 'mean field');
fprintf('%-26s %10s %10s\n', '', 'c^1/3 F^2/3', 'F');

R = zeros(9, 2);
for j = 1:2
    md = modes{j};
    % --- recalibrate this operator ---------------------------------------
    pp = 34.49;  ka = 0.282;
    for outer = 1:8
        ka = fitka(md, pp, NC, TS);
        for it = 1:30
            L = loadOf('baseline', md, ka, pp, NC, TS);
            if L <= 0, break; end
            np = pp * (L / TGT)^0.75;
            if abs(np-pp)/np < 5e-5, pp = np; break; end
            pp = np;
        end
    end
    % --- evaluate ---------------------------------------------------------
    [Lb, Hb, ~ ] = full_('baseline',  md, ka, pp, NC, T16);
    [Lk, ~ , ~ ] = full_('irf8ko',    md, ka, pp, NC, T16);
    [Ld, ~ , ~ ] = full_('depleted',  md, ka, pp, NC, T16);
    [Ln, ~ , L16n] = full_('nocarrier', md, ka, pp, NC, T16);
    [~ , Hd, ~ ] = full_('depleted',  md, ka, pp, NC, T16);

    R(1,j) = pp;
    R(2,j) = ka;
    R(3,j) = Lb;
    R(4,j) = Hb;
    R(5,j) = Lk / max(Lb,eps);
    R(6,j) = Ld / max(Lb,eps);
    R(7,j) = Hd / max(Hb,eps);
    R(8,j) = L16n;
    r6 = Ld / max(Lb,eps);
    if r6 > 0, R(9,j) = log(r6) / log(D.deplete.nuRel); else, R(9,j) = NaN; end
    C.(md) = struct('perPercent',pp,'k_a',ka,'load',Lb,'host',Hb, ...
                    'irf8',R(5,j),'blz',R(6,j),'hostDepl',R(7,j), ...
                    'inert16',R(8,j),'gamma',R(9,j));
end

lbl = {'perPercent (refitted)','k_a (refitted), 1/day', ...
       'control graft load, %','host load at 4 wk, %', ...
       'Irf8-/- load ratio','BLZ945 load ratio', ...
       'host load, 87% depleted','inert carriers, 16 wk %','exponent gamma'};
meas = {'', '', sprintf('measured %.3f, fitted in both',TGT), ...
        'not measured in this preparation', ...
        sprintf('measured %.3f, predicted in both',D.irf8Load.rel), ...
        sprintf('measured %.3f, fitted left / predicted right',D.deplete.rel), ...
        'not measured as an area', 'not measured', ...
        sprintf('%.2f from the two depletion ratios',D.gamma.val)};
for i = 1:9
    fprintf('%-26s %10.4g %10.4g    %s\n', lbl{i}, R(i,1), R(i,2), meas{i});
end

fprintf(['\nThese two columns are NOT an equally calibrated comparison.  m_seed is\n' ...
         'fitted to the BLZ945 ratio under the discrete operator and is not\n' ...
         'refitted under the mean-field one, and the two operators detect a\n' ...
         'deposit differently.  p210_readouts refits the same three targets under\n' ...
         'each rule (Table S3 of the paper).  There the mass-proportional rule\n' ...
         'cannot reach the BLZ945 ratio, but its closest value, 0.206, and its\n' ...
         'Irf8-/- ratio, 0.208, are both within one standard error of the\n' ...
         'measurements, so the depletion data do not select between the rules.\n' ...
         'The rules differ mainly in the host load and the late-time graft load.\n']);
end

% -------------------------------------------------------------------------
function p = mk(variant, md, ka, pp, NC)
p = p210_params(variant);
p.N = NC;  p.k_a = ka;  p.perPercent = pp;  p.readout = md;
% the mean-field operator declares its threshold as an area too, 0.05% of
% it, exactly as the first version of this model did
p.F_thr = 0.05 * pp;
p = p210_derive(p);
end

function L = loadOf(variant, md, ka, pp, NC, TS)
d = p210_post(p210_solve(mk(variant,md,ka,pp,NC), TS));
L = d.loadPct(end);
end

function [L, H, L16] = full_(variant, md, ka, pp, NC, TS)
d = p210_post(p210_solve(mk(variant,md,ka,pp,NC), TS));
i4 = find(TS==28, 1);  if isempty(i4), i4 = numel(TS); end
L = d.loadPct(i4);  H = d.hostPct(i4);  L16 = d.loadPct(end);
end

function ka = fitka(md, pp, NC, TS)
f = @(k) deplOf(md, k, pp, NC, TS) - 0.80;
lo = 1e-3; hi = 50; flo = f(lo); fhi = f(hi);
if flo*fhi > 0, ka = 0.282; return; end
ka = sqrt(lo*hi);
for it = 1:40
    ka = sqrt(lo*hi); fx = f(ka);
    if abs(fx) < 1e-4, return; end
    if fx*flo > 0, lo = ka; flo = fx; else, hi = ka; end
end
end

function dep = deplOf(md, ka, pp, NC, TS)
sol = p210_solve(mk('baseline',md,ka,pp,NC), TS);  G = sol.G;
dep = 1 - sum(G.Vc(G.inHost).*sol.m(end,G.inHost)')/G.Vhost;
end
