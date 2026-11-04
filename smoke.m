% Smoke test: does the solver run and produce sane numbers?
fprintf('MATLAB %s\n', version);
p = p210_params('baseline');
fprintf('k_a = %.5g  k_r = %.5g  ell_m(clean) = %.4f mm\n', ...
        p.k_a, p.k_r, sqrt(p.D_m/p.lam_m));
tic
sol = p210_solve(p, [0 7 14 28 56 112]);
fprintf('solve ok in %.1f s, nt = %d\n', toc, numel(sol.t));
d = p210_post(sol);
fprintf('\n  t(d)   graft%%    host%%   nu_graft(/mm2)  carrFrac\n');
for k = 1:numel(d.t)
    fprintf('%6.0f  %8.4f %8.3f %12.0f %12.3f\n', ...
            d.t(k), d.loadPct(k), d.hostPct(k), d.nuGraft(k), d.carrFrac(k));
end
fprintf('\nmin m = %.4g, min F = %.4g, min nu = %.4g, max nu = %.4g\n', ...
        min(sol.m(:)), min(sol.F(:)), min(sol.nu(:)), max(sol.nu(:)));
fprintf('m in host (end) = %.4f, m at graft centre (end) = %.4f\n', ...
        sol.m(end,end), sol.m(end,1));
fprintf('Vmig = %.4g mm/day, ell_x = %.4f mm, Lambda = %.3f\n', ...
        d.Vmig, d.ell_x, d.Lambda);
fprintf('onset: shallowest %.1f d, deepest %.1f d, fraction detected %.2f\n', ...
        min(d.onsetT), max(d.onsetT), mean(~isnan(d.onsetT)));

fprintf('\n--- no-carrier control ---\n');
pn = p210_params('nocarrier');
sn = p210_solve(pn, [0 28 112]);
dn = p210_post(sn);
fprintf('graft load at 4 wk = %.5f%%, at 16 wk = %.5f%%\n', dn.loadPct(2), dn.loadPct(3));

fprintf('\n--- lesion ---\n');
tic
Lw = p210_lesion('baseline');
fprintf('lesion solve ok in %.1f s\n', toc);
fprintf('nuWin  0h %.0f  6h %.0f  24h %.0f\n', Lw.nuWin(1), Lw.nuWin(7), Lw.nuWin(end));
fprintf('AbWinRel 0h %.3f  6h %.3f  24h %.3f\n', Lw.AbWinRel(1), Lw.AbWinRel(7), Lw.AbWinRel(end));
disp('SMOKE DONE');
