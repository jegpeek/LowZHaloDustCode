; GALEX-SDSS match

datapath = '~/Dropbox/LowZHaloDustData/'
stpath = '~/Documents/StreamDust/'
info = mrdfits(stpath +'gal_info_dr7_v5_2.fits', 1, hdr)

restore, datapath + 'MPAJHU_dr7_match.sav'

galex = mrdfits(stpath + 'ascdr7match_full_dschimin.fits', 1, hdr)
galex = galex[where(galex.mag_fuv ne (-99))]

close_match_radec, galex.ra, galex.dec, info.ra, info.dec, m1, m2, 5.0/3600., 1, miss1

ug = uniq(m2, sort(m2))
m2 = m2[ug]
galex = galex[ug]

mask = fltarr(n_elements(info))

mask[m2] = mask[m2] + 1

mask[whfg] = mask[whfg] + 1

mm2 = mask[m2]
mwhw = mask[whwise]
mwf = mask[whfg]

gmatch = where(mm2 eq 2)
wh3 = where(mask eq 2)

clr = galex[gmatch].mag_fuv - galex[gmatch].mag_nuv - info[wh3].E_BV_SFD*(2.11 - 5.71*info[wh3].E_BV_SFD)

plot, clr, info[wh3].z, psym=3, xra = [-2, 5]


wh = findgen(n_elements(wh3))
ctwh = n_elements(wh3)

zclr = (info[wh3].z)[wh]
clr = clr[wh]

medsm = 500

nct = n_elements(clr)

mwh = (findgen(nct))[medsm/2:nct-medsm/2-1]

zmed = median((clr)[sort(zclr)], medsm)
loadct, 0, /sil 
window, 1
plot, clr, zclr, psym=3, xra=[0, 10], yra=[0, 0.4]
loadct, 13, /sil
oplot, zmed[mwh], (zclr[sort(zclr)])[mwh], color=ps*255, thick=th1
ord = 7
pf = poly_fit((zclr[sort(zclr)])[mwh], zmed[mwh], ord, yfit=yf)
oplot, yf, (zclr)[sort(zclr)], color=200, thick=th2

residc = clr- total(rebin(reform(zclr, ctwh, 1), ctwh, ord+1)^(rebin(reform(findgen(ord+1), 1, ord+1), ctwh, ord+1))*rebin(pf, ctwh, ord+1), 2)
window, 2
plothist, residc, xg, yg, bin=0.01, xra=[-3, 3]

fitg = gaussfit(xg, yg, gf, nterms=4)

print, gf[2]/sqrt(ctwh)*1000.

nuvfuv = replicate({ra:0., dec:0., z:0., color:0.}, ctwh)

nuvfuv.ra = (info[wh3].ra)[wh]
nuvfuv.dec = (info[wh3].dec)[wh]
nuvfuv.z = (info[wh3].z)[wh]
nuvfuv.color = residc/2.11

mwrfits, nuvfuv, datapath+'NUV-FUV.fits'


end