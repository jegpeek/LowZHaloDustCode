
for wise=0, 1 do begin

if wise then begin
	file_stomp = '../STOMP_OUTPUT/MPA-WISE.fit' 
	my_y_tit = textoidl('Color excess g-W1 [mag]')
endif else begin
	 file_stomp = '../STOMP_OUTPUT/MPA-SDSS.fit'
	 my_y_tit = textoidl('Color excess g-r [mag]')
endelse

file_galaxy = '../DATA/fg_MPAJHU.fits'
file1 = '../DATA/g-W1_nod5.fits'
file2 = '../DATA/pg10.fits'

galaxy = MRDFITS(file_galaxy,1)
gW1 = MRDFITS(file1,1)
pg10 = MRDFITS(file2,1)
img = mrdfits('~/Documents/HVCreddening/imaging_dr7_1sig.fits',1, hdr)

a = MRDFITS(file_stomp,1)
a = a(where(a.physical_separation_mpc lt 1))

if wise then fg = gW1 else fg = pg10
	
rollmed, fg.z, fg.color, 0.003, xz, yc
order=7
pf = poly_fit(xz, yc, order, yfit=yf)
nfg = n_elements(fg)
fg.color=fg.color - total(rebin(reform(fg.z, nfg, 1), nfg, order+1)^(rebin(reform(findgen(order+1), 1, order+1), nfg, order+1))*rebin(pf, nfg, order+1), 2)
rollmed, fg.z, fg.color, 0.003, xz, yc

nboot = 1

errs=fltarr(n_elements(a)) + 0.023

allp = fltarr(nboot, 5)
allin = allp

st = systime(/sec)
for i=0, nboot-1 do begin
	inparms = [2d-3, -2, 2, 0, 2]*randomu(seed, 5)
	allin[i, *] = inparms
	loop_bar, i, nboot
	boot = randomu(seed, n_elements(a))*n_elements(a)
	if i eq 0 then boot = findgen(n_elements(a))
	ys = fg[a[boot].master_index].color
	xs = fltarr(3, n_elements(a))
	xs[0, *] = a[boot].physical_separation_mpc*10. ; at 100 kpc
	xs[1, *] = 10^(a[boot].mass_target - 10.77)
	ssfr = a[boot].ssfr_target
	whlt20 = where(ssfr lt -20, ct)
	if ct ne 0 then ssfr(whlt20) = min(ssfr(where(ssfr gt -20)))
	xs[2, *] = 10^(ssfr + 11)
	    	 
	errs=fltarr(n_elements(a)) + 0.023

	faMARS = {x:xs, y:ys, err:errs}
	allp[i, *] = mpfit('truncplfit_marsnossfr', inparms, functargs=faMARS, /quiet)
endfor
if wise then file = 'allpttrnossfr_WISE1_noboot.sav' else file = 'allpttrnossfr_PG101_noboot	.sav'
save, allp, allin, f=file

endfor

end
