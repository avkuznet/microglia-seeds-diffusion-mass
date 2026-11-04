% Diagnostic: where does the free amyloid pool go?
for k2 = [0 0.5 3.0]
  p = p210_params('baseline'); p.N = 200; p.k_a = 0.05; p.k_2 = k2; p.k_r = 0.06;
  sol = p210_solve(p, [0 7 14 28 56 112]);
  G = sol.G; d = p210_post(sol);
  fprintf('\n===== k_2 = %.3g, k_a = %.3g =====\n', k2, p.k_a);
  fprintf(' t(d)  m_far  m_host  m_graft  F_far   F_graft  graft%%  nu_g  b_far\n');
  for j = 1:numel(sol.t)
    mh = sum(G.Vc(G.inHost).*sol.m(j,G.inHost)')/G.Vhost;
    mg = sum(G.Vc(G.inGraft).*sol.m(j,G.inGraft)')/G.Vgraft;
    Fg = sum(G.Vc(G.inGraft).*sol.F(j,G.inGraft)')/G.Vgraft;
    ng = sum(G.Vc(G.inGraft).*sol.nu(j,G.inGraft)')/G.Vgraft;
    fprintf('%5.0f %6.3f %7.3f %8.3f %7.1f %8.2f %7.4f %5.2f %7.2f\n', ...
      sol.t(j), sol.m(j,end), mh, mg, sol.F(j,end), Fg, d.loadPct(j), ng, sol.b(j,end));
  end
end
disp("DIAG DONE");
