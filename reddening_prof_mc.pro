pro reddening_prof_mc, nboot, wise, backcheck=backcheck,use10=use10, angle=angle


datapath = '~/Dropbox/LowZHaloDustData/'

if wise ne 0 then begin
	file_stomp = datapath + 'MPA-WISE'
endif
if wise eq 0 then begin
	file_stomp = datapath + 'MPA-SDSS'
endif

if keyword_set(backcheck) then	file_stomp = file_stomp + '_REVERSE' 
if keyword_set(angle) then file_stomp = file_stomp + '_angspace'
file_stomp = file_stomp + '.fit'

file_galaxy = datapath + 'fg_MPAJHU.fits'
file1 = datapath +'g-W1_nod5.fits'
file2 = datapath + 'pg10.fits'

galaxy = MRDFITS(file_galaxy,1)
gW1 = MRDFITS(file1,1)
pg10 = MRDFITS(file2,1)
;img = mrdfits('~/Documents/HVCreddening/imaging_dr7_1sig.fits',1, hdr)

a = MRDFITS(file_stomp,1)

photo = 0
if photo then begin
	dmf = 1
	restore, datapath + 'photo_magsrads.sav'
	restore, datapath + 'pg10_dr7_match.sav'
	restore, datapath + 'MPAJHU_dr7_match.sav'
	add_tag, a, 'photo', 0., a_photo
	a = a_photo
	a_photo = 0.
	; r band petrosian, for a start. dmag here is in the sense that larger means the obscurer is brighter
	a.photo =  (prads[dmf, whfg])[a.target_index]
	; note -- if I am really doing these matches correctly, why am I getting these -99s? It's not too many, but still. Worrisome for the accuracy of my matches, which may matter in the smallest separations (?)
	a = a[where( ((pmags[dmf, whfg])[a.target_index] ne -9.99) and ((pmags[dmf, whfg])[a.target_index] ne -9999.0) and ((pmags[dmf, m2])[a.master_index] ne -9.99) and ((pmags[dmf, m2])[a.master_index] ne -9999.0))]
	r90 = a.photo                                                    
	r90 = r90(where(r90 gt 0))                                      
	plothist, alog10(r90/60.), xr90, yr90, bin=0.01, /noplot            
endif



;a = a(where(a.physical_separation_mpc lt 1))


tags = tag_names(a)

if keyword_set(angle) then xidx = reform(where(tags eq 'ANGLE')) else xidx = reform(where(tags eq 'PHYSICAL_SEPARATION_MPC'))

; convert angle from arcseconds to 10 arcminutes
if keyword_set(angle) then a.(xidx) = a.(xidx)/600.


if (wise eq 1) then fg = gW1 
if (wise eq 0) then fg = pg10
if (wise eq 3) then fg = mrdfits(datapath +'u-W1_nod5.fits', 1)

if keyword_set(use10) then begin
	restore, '../../HVCreddening/gal_color.sav'
	fg.color = reform(gcs[*, use10])
endif
	
rollmed, fg.z, fg.color, 0.003, xz, yc
order=7
pf = poly_fit(xz, yc, order, yfit=yf)
nfg = n_elements(fg)
fg.color=fg.color - total(rebin(reform(fg.z, nfg, 1), nfg, order+1)^(rebin(reform(findgen(order+1), 1, order+1), nfg, order+1))*rebin(pf, nfg, order+1), 2)
rollmed, fg.z, fg.color, 0.003, xz, yc


zbuf = 0.015
if keyword_set(backcheck) then a = a[where(a.z_target gt (a.z_background + zbuf))] else a = a[where((a.z_target +zbuf) lt a.z_background)] 

fg.color = fg.color-median(fg.color)
normfrac=0.8
fg.color = fg.color - median(fg[a[where(a.(xidx) gt normfrac*max(a.(xidx)))].master_index].color)

errs=fltarr(n_elements(a)) + 0.023

allp = fltarr(nboot, 4)
allin = allp

alimit = 0.500
a = a[where(a.(xidx) lt alimit)]
amin = 0.01
a = a[where(a.(xidx) gt amin)]
if photo then a = a(where(a.photo lt 5))

st = systime(/sec)


for i=0, nboot-1 do begin
	inparms = [1d-2, -2, 1d-2, 1d0]*(0.5 + randomu(seed, 4)*1)
	allin[i, *] = inparms
	loop_bar, i, nboot
	boot = randomu(seed, n_elements(a))*n_elements(a)
	if i eq 0 then boot = findgen(n_elements(a))
	ys = (fg[a[boot].master_index].color);[sort(randomu(seed, n_elements(a)))]
	xs = fltarr(3, n_elements(a))
	xs[0, *] = a[boot].(xidx)*10. ; at 100 kpc
	xs[1, *] = 10^(a[boot].mass_target - 10.77)
	ssfr = a[boot].ssfr_target
	whlt20 = where(ssfr lt -20, ct)
	if ct ne 0 then ssfr(whlt20) = min(ssfr(where(ssfr gt -20)))
	xs[2, *] = 10^(ssfr + 11)
	    	 
	errs=fltarr(n_elements(a)) + 0.023

	faMARS = {x:xs, y:ys, err:errs}
	allp[i, *] = mpfit('plfit_marsnossfr', inparms, functargs=faMARS, /quiet)
endfor

if keyword_set(backcheck) then check = 'bck' else check = 'fwd'
if keyword_set(use10) then check = check + string(use10, f='(I2.2)')
if (wise ne 0) then file = datapath + 'MC_WISE'+ string(wise, f='(I1.1)')+'_' + check+'.sav' else file = datapath + 'MC_SDSS'+check+'.sav'
save, allp, allin, f=file



end
