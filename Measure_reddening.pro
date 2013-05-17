pro Measure_reddening, wise, fit, rc, dofit=dofit, dohist=dohist, ps=ps, backcheck=backcheck, a=a, delmag=delmag, spectromags=spectromags, a0=a0,zbuf = zbuf, use10=use10, galex=galex, angle=angle, dmin=dmin, dmax=dmax, mmin=mmin, mmax=mmax, smin=smin, smax=smax,r90=r90, gw=gw, phmag=phmag, new10=new10

; WISE: 0 is g-r, 1 is g-W1, 2 is g-W2
; FIT: if dofit is set, this is an output of the fit parameters
; DOFIT: if set, run the fit to the data
; DOHIST: if set, run the original histograms
; PS: if set, output an eps file of the histograms
; RC: sdss filter to compare to WISE, if wise != 0. 1 is g-band, 2 is r-band


if keyword_set(backcheck) then ttl = 'Reversed Test' else ttl= 'Reddening Measurement'

resolve_routine,'display_data'

datapath = '~/Dropbox/LowZHaloDustData/'
stpath = '~/Documents/StreamDust/'
if wise ne 0 then begin
	file_stomp = datapath + 'MPA-WISE'
endif
if wise eq 1 then my_y_tit = textoidl('Color excess g-W1 [mag]')
if wise eq 2 then my_y_tit = textoidl('Color excess g-W2 [mag]')
if wise eq 0 then begin
	file_stomp = datapath + 'MPA-SDSS'
	my_y_tit = textoidl('Color excess g-r [mag]')
endif

if keyword_set(gw) then file_stomp = datapath + 'MPA-GALEX-WISE'
if keyword_set(backcheck) then	file_stomp = file_stomp + '_REVERSE' 
if keyword_set(angle) then file_stomp = file_stomp + '_angspace'
file_stomp = file_stomp + '.fit'

file_galaxy = datapath +'fg_MPAJHU.fits'

if rc eq 0 then begin
	filegw1 = datapath +'u-W1_nod5.fits' 
	filegw2 = datapath + 'u-W2_nod5.fits'
endif

if rc eq 1 then begin
	filegw1 = datapath +'g-W1_nod5.fits' 
	filegw2 = datapath + 'g-W2_nod5.fits'
endif
if rc eq 2 then begin
	filegw1 = datapath +'r-W1_nod5.fits' 
	filegw2 = datapath + 'r-W2_nod5.fits'
endif

if keyword_set(gw) then filegw1 = datapath +'NUV-W1.fits' 

file2 = datapath + 'pg10.fits'

if keyword_set(galex) then begin
	file_stomp = datapath + 'MPA-GALEX'
	if keyword_set(backcheck) then	file_stomp = file_stomp + '_REVERSE'
	file_stomp = file_stomp + '.fit'
	fg = mrdfits(datapath +'NUV-FUV.fits')
endif

;if choice eq 0 then begin
galaxy = MRDFITS(file_galaxy,1, /sil)
print, filegw1
if wise eq 1 then fg = MRDFITS(filegw1,1, /sil)
if wise eq 2 then fg = MRDFITS(filegw2,1, /sil)
if wise eq 0 then fg = MRDFITS(file2,1, /sil)

if wise eq 0 then begin
	kc = mrdfits('../../HVCreddening/kcorr0_v4.2_dr7_1sig.fits', 1, hdr, /sil)
	amag = kc.absmag[2]
endif

;if wise ne 0 then begin
;	kc = mrdfits('../../HVCreddening/kcorr0_v4.2_dr7_1sig.fits', 1, hdr, /sil)
;	amag = kc[fg.index].absmag[2]
;endif

if keyword_set(galex) then begin
	fg = mrdfits(datapath +'galex_match_coords.fits', 1)
	whclr = where(fg.color gt 5000, ctclr)
	if ctclr ne 0 then fg[whclr].color = median(fg.color)
endif

;mag =  mrdfits(datapath +'magnitudes_dr7_1sig.fits',1, hdr)
;endif

print, 'file_stomp = ' + file_stomp
if ~keyword_set(a0) then a = MRDFITS(file_stomp,1, /sil)
;a = a0

if keyword_set(use10) then begin
	if keyword_set(new10) then begin
		restore,datapath +'ngcs.sav' 
		fg.color = reform(ngcs[*, use10])
		inf = where(finite(fg.color) eq 0, ctinf)
		if ctinf ne 0 then fg[inf].color = 0
	endif else begin	
		restore, '../../HVCreddening/gal_color.sav'
		fg.color = reform(gcs[*, use10])
	endelse
	c1 = [1,2,0,3,0,1,1,2, 0, 0]
	c2 = [2,3,1,4,2,3,4,4, 3, 4]
	clrs = ['u', 'g', 'r', 'i', 'z']
	fc = [5.15500   ,   3.79300  ,    2.75100   ,   2.08600   ,   1.47900]
	my_y_tit = textoidl('Effective E(B-V): '+clrs[c1[use10]] + '-' + clrs[c2[use10]] +'/' + strcompress(string(fc[c1[use10]]-fc[c2[use10]]),/rem) + ' [mag]')
endif

photo = 0
if photo then begin
	dmf = 1
	restore, datapath + 'photo_magsrads.sav'
	restore, datapath + 'pg10_dr7_match.sav'
	restore, datapath + 'MPAJHU_dr7_match.sav'
	add_tag, a, 'photo', 0., a_photo
	a = a_photo
	a_photo = 0.
	add_tag, a, 'photo2', 0., a_photo
	a = a_photo
	a_photo = 0.
	
	if keyword_set(gw) then begin
		ds = mrdfits(stpath + 'ascdr7match_full_dschimin.fits', 1, hdr)
		
		close_match_radec,ds.ra,ds.dec,galaxy.ra,galaxy.dec,m1t,m2t,4./3600.,1,miss1	
		nuv_mag_target = fltarr(n_elements(galaxy)) + 99
		nuv_mag_target[m2t] = ds[m1t].mag_nuv

	;	close_match_radec,ds.ra,ds.dec,fg.ra,fg.dec,m1b,m2b,4./3600.,1,miss1	
	;	nuv_mag_background = fltarr(n_elements(fg)) + 99
	;	nuv_mag_background[m2b] = ds[m1b].mag_nuv
		
	;	a.photo =  nuv_mag_target[a.target_index] - fg[a.master_index].nuv
		a.photo =  fg[a.master_index].nuv
		a.photo2 = nuv_mag_target[a.target_index]
;		plot, a.photo, (pmags[0, whfg])[a.target_index]
	endif else begin
			; r band petrosian, for a start. dmag here is in the sense that larger means the obscurer is brighter
		a.photo =  (prads[dmf, whfg])[a.target_index]
	
		; note -- if I am really doing these matches correctly, why am I getting these -99s? It's not too many, but still. Worrisome for the accuracy of my matches, which may matter in the smallest separations (?)
		a = a[where( ((pmags[dmf, whfg])[a.target_index] ne -9.99) and ((pmags[dmf, whfg])[a.target_index] ne -9999.0) and ((pmags[dmf, m2])[a.master_index] ne -9.99) and ((pmags[dmf, m2])[a.master_index] ne -9999.0))]
		r90 = a.photo                                                    
		r90 = r90(where(r90 gt 0))                                      
		plothist, alog10(r90/60.), xr90, yr90, bin=0.01, /noplot            
	endelse
endif


; all the tags associated with the stomp output
tags = tag_names(a)

if keyword_set(angle) then xidx = reform(where(tags eq 'ANGLE')) else xidx = reform(where(tags eq 'PHYSICAL_SEPARATION_MPC'))

; convert angle from arcseconds to 10 arcminutes
if keyword_set(angle) then a.(xidx) = a.(xidx)/600.


north=0
minra = 100
maxra = 280

if keyword_set(spectromags) then begin

	restore, datapath + 'pg10_spectrormi.sav', /ver
	print, 'SMS'
	fg.color = spectro_clr
endif

docorrect=1
; option to do a polynomial correction in z. SUGGESTED FOR CURRENT REDUCTION
if docorrect then begin
	
	rollmed, fg.z, fg.color, 0.003, xz, yc
	order=7
	pf = poly_fit(xz, yc, order, yfit=yf)
	nfg = n_elements(fg)
	fg.color=fg.color - total(rebin(reform(fg.z, nfg, 1), nfg, order+1)^(rebin(reform(findgen(order+1), 1, order+1), nfg, order+1))*rebin(pf, nfg, order+1), 2)
	rollmed, fg.z, fg.color, 0.003, xz, yc
;	plot, xz, yc, psym=4, color=100	
endif

if keyword_set(dohist) then begin	
	if north then a = a[ where(fg[a.master_index].ra gt minra and fg[a.master_index].ra lt maxra)]
	; set parameters for limits on mass, SSFR. 14, 1, -20, -1 is effectively without limits
	if keyword_set(backcheck) then a = a[where(a.z_target gt (a.z_background + zbuf))] else a = a[where((a.z_target +zbuf) lt a.z_background)] 


	if photo and not keyword_Set(gw) then begin	
		r90 = a.photo                                                    
		r90 = r90(where(r90 gt 0))                                      
		plothist, alog10(r90/60.), xr90, yr90, bin=0.01, /noplot            
	endif

	if ~keyword_set(mmax) then	mmax = 14
	if ~keyword_set(mmin) then	mmin = 1
	if ~keyword_set(smin) then	smin = (-20.0)
	if ~keyword_set(smax) then	smax = (-1.0)
	if ~keyword_set(dmin) then dmin =0		
	if ~keyword_set(dmax) then dmax = 10000 
	a = a[where(a.mass_target lt mmax and a.mass_target gt mmin)]
	a = a[where(a.ssfr_target lt smax and a.ssfr_target gt smin)]
	print, 'median ssfr target in STOMP = '+ string(median(a.ssfr_target))
	if photo then a = a[where(a.photo lt dmax and a.photo gt dmin)]
	if photo then print, median(a.photo)
	print, 'median mass target in STOMP = ' + string(median(a.mass_target))
;	if wise ne 0 then delmag = a.pmagr_target - fg[a.master_index].dered_mag[2]
;	if wise eq 0 then delmag = a.pmagr_target - mag[a.master_index].dered_mag[2]
;	whdm = where(delmag gt min(delmag))
; 	a = a[whdm]
; 	delmag = delmag[whdm]
    n_r_bin = 10
    r_vector = make_vector(0.02,3. < (max(a.(xidx))), /log,n_r_bin)
    color_single = {mean:0.,mean_err:0.,count:0L,median:0., medbterr:0., meanbterr:0., medzbin:0.}
    color_list = REPLICATE(color_single, n_r_bin)
    
	mdclr=  median(fg.color)
	mnclr=  mean(fg.color)
	;hgz = h2d_ri(fg.z, amag, zvec=fg.color-mnclr, 0.02, 0.2, xrng=[0, 0.3], yrng=[-24, -16], zimg=clrh2d)
	
	if north then begin
		mdclr=  median(fg[where(fg.ra gt minRA and fg.ra lt maxRA)].color)
		mnclr=  mean(fg[where(fg.ra gt minRA and fg.ra lt maxRA)].color)
	endif
	donormhist = 0
    for i_bin=0,n_r_bin-1 do begin
    	loop_bar, i_bin, n_r_bin
        ind_in_bin = where( ((a.(xidx)) gt r_vector[i_bin].bound_min) AND $
                            ((a.(xidx)) lt r_vector[i_bin].bound_max), ct)
		if donormhist then begin
			hgz_ind_in_bin = h2d_ri(a[ind_in_bin].z_background, amag[a[ind_in_bin].master_index], 0.02, 0.2, xrng=[0, 0.3], yrng=[-24, -16])
			
			;hgind = histogram(fg[a[ind_in_bin].master_index].z, min=0.0, max=0.3-1d-6, bin=0.01)
			;clrind = hgind*0.0
			;for j=0, n_elements(hgind)-1 do begin
			;	if hgind[j] ne 0 then clrind[j] = median(fg[(ri[ri[j]:ri[j+1]-1])[randomu(seed, hgind[j]*100)*hgz[j]]].color)-mdclr				
			;endfor
			color_list[i_bin].medzbin = total(hgz_ind_in_bin*clrh2d)/total(hgz_ind_in_bin)
		endif
		
        color_list[i_bin].count = n_elements(ind_in_bin)
        color_list[i_bin].mean = mean(fg[a[ind_in_bin].master_index].color)-mnclr
        if ct ne 1 then color_list[i_bin].median = MEDIAN(fg[a[ind_in_bin].master_index].color)-mdclr else color_list[i_bin].median = fg[a[ind_in_bin].master_index].color-mdclr
        color_list[i_bin].mean_err = STDDEV(fg[a[ind_in_bin].master_index].color)/sqrt(n_elements(ind_in_bin))
        
		nb = 0
		if nb gt 0 then begin
			medvals = fltarr(nb)
			meanvals = fltarr(nb)
			for j=0, nb-1 do begin
				rnd = randomu(seed, ct)*ct
				medvals[j] = MEDIAN(a[ind_in_bin[rnd]].color)
				meanvals[j] =MEAN(a[ind_in_bin[rnd]].color)
			endfor        
			color_list[i_bin].meanbterr = STDDEV(meanvals)
			color_list[i_bin].medbterr = STDDEV(medvals)
		endif
    endfor
	
	circle, /fill
    x = r_vector.mean_2d*1000
    print, x
    inorm = color_list[n_r_bin-1].median
    ; subtracting off some kind of error from the redshift distribution?
    y = color_list.median-inorm
    print, 'lastbin='
    print, color_list[n_r_bin-1].median
    y_err = color_list.mean_err
    print,y
    print,color_list.count
	th = 5
    my_Yr = [1e-5,1e-1]
    my_xr = [0.02,3.0]*1000.
	cc = ['g-r', 'g-W1', 'g-W2']
    if keyword_set(ps) then psopen, datapath+'reddening_'+cc[wise]+'m' + string(mmin, f='(F4.1)') +'--'+ string(mmax, f='(F4.1)') + 's' + string(smin*(-1), f='(F4.1)') +'--'+ string(smax*(-1), f='(F4.1)'), /helvetica, xsi=4.5, ysi=3, /inches, /color, /encapsulated
	if ~keyword_set(ps) then ps=0
    loadct, 0
    if ps then !p.font=0 else !p.font = (-1)
	if keyword_set(angle) then begin
		xtit='separation [arcmin]'
    	xscl = 1d-2
    endif else begin
    	xtit='separation [kpc]'
    	xscl=1
    endelse
    plot,x*xscl,y,/xlog,psym=4,yr=my_yr,ylog=1,xr=my_xr*xscl,$
      xtit=xtit,$
      ytit=my_y_tit, /xs, thick=th, xthick=th, ythick=th, /nodata, title=ttl, syms=0.6
      
   	if photo and (not keyword_set(gw)) then begin
		r90 = a.photo                                                    
		r90 = r90(where(r90 gt 0))                                      
		
		whang = where((10^xr90 gt dmin/60.) and (10^xr90 lt dmax/60.), ctang)
		polyfill, [(10^xr90)[whang[0]], (10^xr90)[whang] > 10^!x.crange[0], (10^xr90)[whang[ctang-1]]] ,[10^!y.crange[0], (yr90*0.1/max(yr90))[whang],10^!y.crange[0]] > 10^!y.crange[0], color=200
		oplot, 10^xr90, yr90*0.1/max(yr90), psym=10, thick=2 
		xyouts, 0.1, 0.8, 'Petro R90 histogram', /normal
	endif
	
;    my_oploterr,x,y,y_err,psym=4,miny=my_yr[0], thick=th
    print,y
    oplot,x*xscl,y,psym=8, color=getcolor('red',1), syms=0.6
        my_oploterr,x*xscl,y,y_err,psym=4,miny=my_yr[0],errcolor=getcolor('red',1), thick=th
    ylast = y[n_elements(y)-1]
    ;oplot, x, color_list.medzbin,color=getcolor('red',1), psym=-2, line=1	
    ;oplot, x, color_list.medzbin*(-1),color=getcolor('blue',1), psym=-2
	print,color_list.medzbin
	oplot,x*1.05*xscl,y*(-1),psym=8, color=getcolor('blue',1), syms=0.6     
	my_oploterr,x*1.05*xscl,y*(-1),y_err,psym=4,miny=my_yr[0],errcolor=getcolor('blue',1), linestyle=1, thick=th
 
  ;  oplot,x*1.05,y-ylast,psym=4, color=getcolor('green',1)
   ; 	my_oploterr,x*1.05,y-ylast,y_err,psym=4,miny=my_yr[0],errcolor=getcolor('green',1)
	
	xax = alog10(findgen(100)*10)
	;oplot, xscl*10^xax, 0.5*4.14d-3/3.1*((10.^xax)/100.)^(-0.84), color=200, thick=2*th, lines=1
	;oplot, xscl*10^xax, 4.14d-3/3.1*((10.^xax)/100.)^(-0.84), color=200, thick=2*th, lines=2
	;oplot, xscl*10^xax, 5*4.14d-3/3.1*((10.^xax)/100.)^(-0.84), color=200, thick=2*th, lines=1

	;xyouts, 0.7, 0.8, mean(10^(a.mass_target)), /norm, charsize=3-ps*2
	;xyouts, 0.7, 0.7, mean(a.z_target), /norm, charsize=3-ps*2
	print, 	(mean(10^(a.mass_target)))
	if keyword_set(ps) then psclose
endif


; do the actual fit to the data
if keyword_set(dofit) then begin
	
	if photo then a = a[where(a.photo lt dmax and a.photo gt dmin)]
	
	
	if keyword_set(backcheck) then a = a[where(a.z_target gt (a.z_background + zbuf))] else a = a[where((a.z_target +zbuf) lt a.z_background)] 

	; option to use median fitting. SUGGESTED.
	usemed = 1
	; brice's fix to deal with flattening population with z.
	impactnorm=1
	fg.color = fg.color-median(fg.color)
	if keyword_set(dohist) then fg.color = fg.color - inorm else begin
		normfracmin = 0.9/max(a.(xidx))
		normfracmax=1/max(a.(xidx))
		rmax = 0.848
		rmin = 0.573
		if usemed then begin
			method = 'plfit_mars'
			fg.color = fg.color-median(fg.color)
			print, 'normfrac median'
			print, median(fg[a[where( (a.(xidx) gt rmin) and (a.(xidx) lt rmax))].master_index].color)
			if impactnorm then begin
				inorm = median(fg[a[where( (a.(xidx) gt rmin) and (a.(xidx) lt rmax))].master_index].color)
				print, 'inorm = ' + string(inorm)
				fg.color = fg.color - inorm
			endif
		endif else begin
			method = 'plfit'
			fg.color = fg.color-mean(fg.color)
			if impactnorm then fg.color = fg.color - mean(fg[a[where( (a.(xidx) gt rmin) and (a.(xidx) lt rmax))].master_index].color)
		endelse
	endelse
	
;	print, max(a.(xidx))
	alimit = 1.0
	a = a[where(a.(xidx) lt alimit)]
	amin = 0.02
	a = a[where(a.(xidx) gt amin)]
	
	if usemed then method = 'plfit_mars' else 	method = 'plfit'
	
	
	; option to avoid fitting to specific star formation rate	
	nossfr=1
	if nossfr then method = method+'nossfr'

	trunc = 0
	
	if keyword_set(gw) and keyword_set(angle) then ys = fg[a.master_index].color - compute_fbb(a.photo, a.angle*10.) else ys = fg[a.master_index].color
	xs = fltarr(3, n_elements(a))
	xs[0, *] = a.(xidx)*10. ; in units of 100 kpc
	xs[1, *] = 10^(a.mass_target - 10.77) ; in units of 6 x 10^10 solar masses
	ssfr = a.ssfr_target
	czsfr = where(ssfr lt -20, ct)
	if ct ne 0 then ssfr(czsfr) = min(ssfr(where(ssfr gt -20))) ; get rid of a few crazy outliers
	xs[2, *] = 10^(a.ssfr_target + 11) ; in units of the ~median SSFR, 10^-11 
	; HACK for crayon brightness (?)
	
	errs=fltarr(n_elements(a)) + 0.023 ; this doesn't actually get used in the median fit, FYI
	faMARS = {x:xs, y:ys, err:errs}
		
	; truncation with a single radius (1) or mass dependency (2)
	if trunc eq 1 then method = 'trunc' + method
	if trunc eq 2 then method = 'trunc2' + method
	inparms = [-1d-2, -1d0, 1d0, 1d0];*randomu(seed, 4)
	if ~nossfr then inparms = [inparms, 1]
	if trunc eq 1 then inparms = [inparms , 1]
	if trunc eq 2 then inparms = [inparms , 1, 0.1]
	print, 'method = ' + method
	fit = mpfit(method, inparms, functargs=faMARS)

	if photo then begin
		print, 'median mag:'
		print, median(a.photo)
		phmag = median(a.photo)
	endif
	rollmed, alog10(a.(xidx)*10.), ys, 0.192, xc, yc
	oplot, 10^xc, (-1)*yc, psym=-1, color=getcolor('blue',1)
	oplot, 10^xc, yc, psym=-1, color=getcolor('red',1)
	xax = findgen(10000)/3.
	oplot, xax, fit[0]*(-1)*((xax/100)^fit[1]), color=getcolor('blue',1)
	oplot, xax, fit[0]*((xax/100)^fit[1]), color=getcolor('red',1)

	stop
	if photo then rollmed, alog10(a.(xidx)*10.), a.photo2, 0.192, xc, yc

endif


end
