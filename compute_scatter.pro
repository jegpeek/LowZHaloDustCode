pro compute_scatter, allp, R_fid, mass_lstar, ssfr_11, value, err_rng

p0 = reform(allp)
sz = size(p0)
if sz[2] eq 4 then p0 = [[p0], [fltarr(300)+1e6]]
y = p0[*, 0]*(R_fid/100.)^(p0[*, 1])*mass_lstar^(p0[*, 2])*ssfr_11^(p0[*, 3])*exp((-1)*(R_fid/100.)/p0[*, 4])

value = y[0]

boots = y[1:sz[1]-1]

boots = boots[sort(boots)]

onesig_min = 0.5 - 0.682689492137/2.
onesig_pls = 0.5 + 0.682689492137/2.

err_rng=fltarr(2)
err_rng[0] = boots[onesig_min*(sz[1]-1)]
err_rng[1] = boots[onesig_pls*(sz[1]-1)]

end 