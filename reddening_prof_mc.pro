pro reddening_prof_mc, nboot, backcheck=backcheck,use10=use10


wise=0
datapath = '~/Dropbox/LowZHaloDustData/'
if wise ne 0 then begin
	file_stomp = datapath + 'MPA-WISE.fit'
	if keyword_set(backcheck) then file_stomp = datapath + 'MPA-WISE_REVERSE.fit' 
endif
if wise eq 0 then begin
	 file_stomp = datapath + 'MPA-SDSS.fit'
	 if keyword_set(backcheck) then begin
	 	file_stomp = datapath + 'MPA-SDSS_REVERSE.fit' 
	endif
	 my_y_tit = textoidl('Color excess g-r [mag]')
endif

file_galaxy = datapath + 'fg_MPAJHU.fits'
file1 = datapath +'g-W1_nod5.fits'
file2 = datapath + 'pg10.fits'



galaxy = MRDFITS(file_galaxy,1)
gW1 = MRDFITS(file1,1)
pg10 = MRDFITS(file2,1)
;img = mrdfits('~/Documents/HVCreddening/imaging_dr7_1sig.fits',1, hdr)

a = MRDFITS(file_stomp,1)
a = a(where(a.physical_separation_mpc lt 1))

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


zbuf = 0.02
if keyword_set(backcheck) then a = a[where(a.z_target gt (a.z_background + zbuf))] else a = a[where((a.z_target +zbuf) lt a.z_background)] 

fg.color = fg.color-median(fg.color)
normfrac=0.8
fg.color = fg.color - median(fg[a[where(a.physical_separation_mpc gt normfrac*max(a.physical_separation_mpc))].master_index].color)

errs=fltarr(n_elements(a)) + 0.023

allp = fltarr(nboot, 4)
allin = allp


	alimit = 0.500
	a = a[where(a.physical_separation_mpc lt alimit)]
	amin = 0.01	
	a = a[where(a.physical_separation_mpc gt amin)]


st = systime(/sec)


for i=0, nboot-1 do begin
	inparms = [1d-2, -2, 1d-2, 1d0]
	allin[i, *] = inparms
	loop_bar, i, nboot
	boot = randomu(seed, n_elements(a))*n_elements(a)
	if i eq 0 then boot = findgen(n_elements(a))
	ys = (fg[a[boot].master_index].color);[sort(randomu(seed, n_elements(a)))]
	xs = fltarr(3, n_elements(a))
	xs[0, *] = a[boot].physical_separation_mpc*10. ; at 100 kpc
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
if (wise ne 0) then file = 'MC_WISE'+ string(wise, f='(I1.1)')+'_' + check+'.sav' else file = 'MC_SDSS'+check+'.sav'
save, allp, allin, f=file



end
