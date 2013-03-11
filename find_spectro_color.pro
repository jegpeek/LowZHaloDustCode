; what is the scatter in the spectro-derived colors in g-r?

;gi = mrdfits('gal_info_dr7_v5_2.fits', 1, hdr)
;restore, 'pg10togal_info.sav'

rollmed, gi[indx].z, gi[indx].spectro_mag[0] -  gi[indx].spectro_mag[1],  0.01, xs, ys

order=7
pf = poly_fit(xs[1:28], ys[1:28], order, yfit=yf)
nfg = n_elements(indx)
gmr = gi[indx].spectro_mag[0] -  gi[indx].spectro_mag[1]
gmr_z = gmr - total(rebin(reform(gi[indx].z, nfg, 1), nfg, order+1)^(rebin(reform(findgen(order+1), 1, order+1), nfg, order+1))*rebin(pf, nfg, order+1), 2)

plot, gi[indx].z, gmr_z, psym=3, yra=[-0.2, 0.2]

plothist, gmr_z, bin=0.001, xx, yy, xra=[-1, 1]

fitg = gaussfit(xx[0:1999], yy[0:1999], gf, nterms=4)

amagr_spectro = gi[indx].spectro_mag[1] - 5*(alog10(lumdist(gi[indx].z)*1d6) - 1)

rollmed, amagr_spectro,gmr_z,  0.1, xs, ys
whsub = where(xs gt (-22.9) and xs lt (-16.8))
order=7
pf = poly_fit(xs[whsub], ys[whsub], order, yfit=yf)
nfg = n_elements(indx)
gmr_za = gmr_z - total(rebin(reform(amagr_spectro, nfg, 1), nfg, order+1)^(rebin(reform(findgen(order+1), 1, order+1), nfg, order+1))*rebin(pf, nfg, order+1), 2)

plothist, gmr_za, bin=0.001, xx, yy, xra=[-1, 1]

fitg = gaussfit(xx[0:1999], yy[0:1999], gf, nterms=4)

rollmed, gi[indx].z, gmr_za,  0.01, xs, ys

save, gmr_z, f='pg10_spectrogmr.sav'

end

