; GALEX-SDSS-WISE match

datapath = '~/Dropbox/LowZHaloDustData/'
stpath = '~/Documents/StreamDust/'
info = mrdfits(stpath +'gal_info_dr7_v5_2.fits', 1, hdr)

restore, datapath + 'MPAJHU_dr7_match.sav'

galex = mrdfits(stpath + 'ascdr7match_full_dschimin.fits', 1, hdr)
;galex = galex(where(galex.mag_fuv ne (-99)))
close_match_radec, galex.ra, galex.dec, info.ra, info.dec, m1, m2, 5.0/3600., 1, miss1

dr7w = mrdfits(stpath + 'dr7wise.fits', 1, hdr)  

whwise = dr7w.cntr_u-1l

; uniquify
uw = uniq(whwise, sort(whwise))
whwise = whwise[uw]
dr7w = dr7w[uw]

; uniquify
ug = uniq(m2, sort(m2))
m2 = m2[ug]
galex = galex[ug]

dooiiha = 0
 begin
	oiiha = mrdfits('~/Documents/HVCreddening/oII_halpha.fits', 1)
	hash_oiiha = (long(oiiha.plate) +long(oiiha.fiber)*1000l)*oiiha.mjd
	hash_whfg = (long(info[whfg].plateid) +long(info[whfg].fiberid)*1000l)*info[whfg].mjd
	match2, hash_whfg, hash_oiiha, suba, subb
endif

mask = fltarr(n_elements(info))

; is wise
mask[whwise] = mask[whwise] + 1
; is non duplicated

if dooiiha then mask[whfg[subb]] = mask[whfg[subb]] + 1 else mask[whfg] = mask[whfg] + 1
; is galex
mask[m2] = mask[m2] + 1
; has oiiha
;mask[m2oh] = mask[m2oh] +1


mm2 = mask[m2]
mwhw = mask[whwise]
mwf = mask[whfg]

gmatch = where(mm2 eq 3)
wmatch = where(mwhw eq 3)
smatch = where(mwf eq 3)

wh3 = where(mask eq 3)

clr = galex[gmatch].mag_nuv - dr7w[wmatch].W1MPRO - info[wh3].E_BV_SFD*9
loadct, 0, /sil

; no liners or non starformers:
;wh = where(((oiiha[suba[smatch]].oii_ew gt 5) or (oiiha[suba[smatch]].halpha_ew gt 2)) and (oiiha[suba[smatch]].halpha_ew gt 0.5*oiiha[suba[smatch]].oii_ew), ctwh)

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
; seems like the OII matching is hardly worth it -- will do later, but for now, just do the full match

nuvw1 = replicate({ra:0., dec:0., z:0., color:0., nuv:0.}, ctwh)

nuvw1.ra = (info[wh3].ra)[wh]
nuvw1.dec = (info[wh3].dec)[wh]
nuvw1.z = (info[wh3].z)[wh]
nuvw1.color = residc/8.
nuvw1.nuv = galex[gmatch].mag_nuv
stop
mwrfits, nuvw1, datapath + 'NUV-W1.fits'

end