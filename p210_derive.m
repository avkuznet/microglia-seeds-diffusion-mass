function p = p210_derive(p)
%P210_DERIVE  Fill in the quantities that follow from the declared ones.
%
%   Called at the end of p210_params, and again by p210_calibrate whenever
%   it changes p.perPercent, so that the two never disagree.
%
%   kap is a property of the stain and the microscope and is fixed once,
%   from the reference host: p.F_hRef_pct of area carried by p.c_hRef
%   deposits per mm^3 holding p.F_hRef_pct * p.perPercent of mass.  It must
%   not be re-derived per variant, or each variant would be read with a
%   different microscope.
%
%   Tissue declared to carry a different area - the aged host, the counting
%   window of the lesion experiment - is taken to carry the same deposits in
%   different numbers rather than deposits of a different size, so both the
%   number and the mass scale with the declared area.  The readout is then
%   exactly linear in the declared area, as it has to be.

%   The reference host is read through the same visibility factor as
%   everything else, so that factor enters kap too; otherwise the host would
%   not read back the area it was declared to have.
p.F_hRef = p.F_hRef_pct * p.perPercent;
visRef   = 1 / (1 + (p.detRatio * p.m_seed / (p.F_hRef / p.c_hRef))^3);
p.kap    = p.F_hRef_pct / (p.c_hRef^(1/3) * p.F_hRef^(2/3) * visRef);

rh      = p.F_h_pct / p.F_hRef_pct;
p.F_h   = p.F_hRef * rh;
p.c_h   = p.c_hRef * rh;
p.c_sat = p.c_h;
p.m_c   = p.F_c / p.c_h;

if isfield(p,'F_hole_pct')
    rw       = p.F_hole_pct / p.F_h_pct;
    p.F_hole = p.F_h * rw;
    p.c_hole = p.c_h * rw;
end
end
