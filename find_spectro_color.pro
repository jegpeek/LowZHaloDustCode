; what is the scatter in the spectro-derived colors in g-r?

gi = mrdfits('../gal_info_dr7_v5_2.fits', 1, hdr)
restore, 'pg10togal_info.sav'

fc_gri = [3.79300  ,    2.75100   ,   2.08600]


clr = ['g', 'r', 'i']
bluefilter = 0
redfilter = 2

rollmed, gi[indx].z, gi[indx].spectro_mag[bluefilter] -  gi[indx].spectro_mag[redfilter],  0.01, xs, ys

order=7
pf = poly_fit(xs[1:28], ys[1:28], order, yfit=yf)
nfg = n_elements(indx)
gmr = gi[indx].spectro_mag[bluefilter] -  gi[indx].spectro_mag[redfilter]
gmr_z = gmr - total(rebin(reform(gi[indx].z, nfg, 1), nfg, order+1)^(rebin(reform(findgen(order+1), 1, order+1), nfg, order+1))*rebin(pf, nfg, order+1), 2)

plot, gi[indx].z, gmr_z, psym=3, yra=[-0.2, 0.2]

plothist, gmr_z, bin=0.001, xx, yy, xra=[-1, 1]

fitg = gaussfit(xx[0:1999], yy[0:1999], gf_z, nterms=4)

print, gf_z[2]/(fc_gri[bluefilter]-fc_gri[redfilter])

amagr_spectro = gi[indx].spectro_mag[redfilter] - 5*(alog10(lumdist(gi[indx].z)*1d6) - 1)

rollmed, amagr_spectro,gmr_z,  0.1, xs, ys
whsub = where(xs gt (-22.9) and xs lt (-16.8))
order=7
pf = poly_fit(xs[whsub], ys[whsub], order, yfit=yf)
nfg = n_elements(indx)
gmr_za = gmr_z - total(rebin(reform(amagr_spectro, nfg, 1), nfg, order+1)^(rebin(reform(findgen(order+1), 1, order+1), nfg, order+1))*rebin(pf, nfg, order+1), 2)

plothist, gmr_za, bin=0.001, xx, yy, xra=[-1, 1]

fitg = gaussfit(xx[0:1999], yy[0:1999], gf_za, nterms=4)

rollmed, gi[indx].z, gmr_za,  0.01, xs, ys

spectro_clr = gmr_z/(fc_gri[bluefilter]-fc_gri[redfilter])

save, spectro_clr, f='pg10_spectro'+clr[bluefilter]+'m'+clr[redfilter]+'.sav'

end

