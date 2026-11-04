function C = p210_competition()
%P210_COMPETITION  The two observation operators, side by side.
%
%   C = P210_COMPETITION
%
%   The histological readout is the area that amyloid deposits cover.  The
%   first version of this model took that area to be proportional to the
%   deposit MASS in a voxel, gated by a threshold on it, which is a mean
%   field over a deposit population whose number never enters.  Reviewers of
%   that version asked for sub-voxel competition between deposits, and they
%   were right to, but the place it is needed turns out to be the readout
%   rather than the kinetics.
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
%   This routine recalibrates each operator on its own terms - k_a to the
%   measured host monomer depletion, perPercent to the measured 4-week graft
%   load - so that the comparison is like for like, and then reports what
%   each one predicts for the three perturbations that were not fitted.

D  = p210_data();
TS = [0 7 14 21 28];
T16 = [0 7 14 28 112];
NC = 240;
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
    for outer = 1:4
        ka = fitka(md, pp, NC, TS);
        for it = 1:10
            L = loadOf('baseline', md, ka, pp, NC, TS);
            if L <= 0, break; end
            np = pp * (L / TGT)^0.75;
            if abs(np-pp)/np < 1e-3, pp = np; break; end
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
meas = {'', '', sprintf('%.3f, fitted in both',TGT), '~2 (Grathwohl 2009)', ...
        sprintf('%.3f, predicted in both',D.irf8Load.rel), ...
        sprintf('%.3f, fitted left / predicted right',D.deplete.rel), ...
        '1.07 (ED Fig. 6b)', '0', sprintf('%.2f',D.gamma.val)};
for i = 1:9
    fprintf('%-26s %10.4g %10.4g    measured %s\n', lbl{i}, R(i,1), R(i,2), meas{i});
end

fprintf(['\nThe deposit number is what the depleted graft is short of.  Leaving\n' ...
         'it out of the readout makes the load respond to the carrier count far\n' ...
         'more steeply than it does, because the mass is then spread thin over a\n' ...
         'deposit population the operator assumes is still there, and the\n' ...
         'threshold cuts it twice.  Putting it in gives an exponent close to the\n' ...
         'measured one and leaves the seeding bottleneck intact, now for a\n' ...
         'physical reason rather than an assumed threshold: with inert carriers\n' ...
         'there are no centres for the amyloid to collect on, so what little\n' ...
         'forms is spread over specks too small to resolve.\n']);
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
for it = 1:22
    ka = sqrt(lo*hi); fx = f(ka);
    if abs(fx) < 3e-3, return; end
    if fx*flo > 0, lo = ka; flo = fx; else, hi = ka; end
end
end

function dep = deplOf(md, ka, pp, NC, TS)
sol = p210_solve(mk('baseline',md,ka,pp,NC), TS);  G = sol.G;
dep = 1 - sum(G.Vc(G.inHost).*sol.m(end,G.inHost)')/G.Vhost;
end
