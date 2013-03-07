pro Measure_reddening, wise, fit, dofit=dofit, dohist=dohist, ps=ps

; WISE: 0 is g-r, 1 is g-W1, 2 is g-W2
; FIT: if dofit is set, this is an output of the fit parameters
; DOFIT: if set, run the fit to the data
; DOHIST: if set, run the original histograms
; PS: if set, output an eps file of the histograms

resolve_routine,'display_data'

datapath = '~/Dropbox/LowZHaloDustData/'
if wise eq 1 then begin
	file_stomp = datapath + 'MPA-WISE.fit' 
	my_y_tit = textoidl('Color excess g-W1 [mag]')
endif
if wise eq 2 then begin
	file_stomp = datapath + 'MPA-WISE.fit' 
	my_y_tit = textoidl('Color excess g-W2 [mag]')
endif
if wise eq 0 then begin
	 file_stomp = datapath + 'MPA-SDSS.fit'
	 my_y_tit = textoidl('Color excess g-r [mag]')
endif

file_galaxy = datapath +'fg_MPAJHU.fits'
filegw1 = datapath +'g-W1_nod5.fits' 
filegw2 = datapath + 'g-W2_nod5.fits'
file2 = datapath + 'pg10.fits'

;if choice eq 0 then begin
galaxy = MRDFITS(file_galaxy,1)

if wise eq 1 then fg = MRDFITS(filegw1,1)
if wise eq 2 then fg = MRDFITS(filegw2,1)
if wise eq 0 then fg = MRDFITS(file2,1)

; don't need this, I don't think
;img =  mrdfits(datapath +'imaging_dr7_1sig.fits',1, hdr)
;endif

trunc = 1
a = MRDFITS(file_stomp,1)
;a = a[where(a.physical_separation_mpc lt 1)]
; option to get rid of the southern strips. NOT SUGGESTED.

north=0
minra = 100
maxra = 280

docorrect=1
; option to do a polynomial correction in z. SUGGESTED FOR CURRENT REDECTION
if docorrect then begin
	rollmed, fg.z, fg.color, 0.003, xz, yc
	order=7
	pf = poly_fit(xz, yc, order, yfit=yf)
	nfg = n_elements(fg)
	fg.color=fg.color - total(rebin(reform(fg.z, nfg, 1), nfg, order+1)^(rebin(reform(findgen(order+1), 1, order+1), nfg, order+1))*rebin(pf, nfg, order+1), 2)
	rollmed, fg.z, fg.color, 0.003, xz, yc
	;oplot, xz, yc, psym=4, color=100
endif


; do the actual fit to the data
if keyword_set(dofit) then begin
; option to use median fitting. SUGGESTED.
	usemed = 1
	if usemed then begin
		method = 'plfit_mars'
		fg.color = fg.color-median(fg.color)
	endif else begin
		method = 'plfit'
		fg.color = fg.color-mean(fg.color)
	endelse
	
	; option to avoid fitting to specific starformation rate	
	nossfr=1
	if nossfr then method = method+'nossfr'
		

	ys = fg[a.master_index].color
	xs = fltarr(3, n_elements(a))
	xs[0, *] = a.physical_separation_mpc*10. ; in units of 100 kpc
	xs[1, *] = 10^(a.mass_target - 10.77) ; in units of 6 x 10^10 solar masses
	ssfr = a.ssfr_target
	ssfr(where(ssfr lt -20)) = min(ssfr(where(ssfr gt -20))) ; get rid of a few crazy outliers
	xs[2, *] = 10^(a.ssfr_target + 11) ; in units of the ~median SSFR, 10^-11 
	errs=fltarr(n_elements(a)) + 0.023 ; this doesn't actually get used in the median fit, FYI
	faMARS = {x:xs, y:ys, err:errs}
		
	; truncation with a single radius (1) or mass dependency (2)
	if trunc eq 1 then method = 'trunc' + method
	if trunc eq 2 then method = 'trunc2' + method
	inparms = [1d-3, -1, 1, 1]
	if ~nossfr then inparms = [inparms, 1]
	if trunc eq 1 then inparms = [inparms , 1]
	if trunc eq 2 then inparms = [inparms , 1, 0.1]
	fit = mpfit(method, inparms, functargs=faMARS)
endif

if keyword_set(dohist) then begin	
	if north then a = a[ where(fg[a.master_index].ra gt minra and fg[a.master_index].ra lt maxra)]
	; set parameters for limits on mass, SSFR. 14, 1, -20, -1 is effectively without limits
	mmax = 14.0
	mmin = 1.0
	smin = (-20.0)
	smax = (-1.0)
	a = a[where(a.mass_target lt mmax and a.mass_target gt mmin)]
	print, median(a.ssfr_target)
	a = a[where(a.ssfr_target lt smax and a.ssfr_target gt smin)]
	zmin = 0.04
	a = a[where(a.z_target gt zmin)]
    n_r_bin = 10
    r_vector = make_vector(0.02,3.,/log,n_r_bin)
    color_single = {mean:0.,mean_err:0.,count:0L,median:0., medbterr:0., meanbterr:0.}
    color_list = REPLICATE(color_single, n_r_bin)
    
	mdclr=  median(fg.color)
	mnclr=  mean(fg.color)
	
	if north then begin
		mdclr=  median(fg[where(fg.ra gt minRA and fg.ra lt maxRA)].color)
		mnclr=  mean(fg[where(fg.ra gt minRA and fg.ra lt maxRA)].color)
	endif
	
    for i_bin=0,n_r_bin-1 do begin
    	loop_bar, i_bin, n_r_bin
        ind_in_bin = where( ((a.physical_separation_mpc) gt r_vector[i_bin].bound_min) AND $
                            ((a.physical_separation_mpc) lt r_vector[i_bin].bound_max), ct)

        color_list[i_bin].count = n_elements(ind_in_bin)
        color_list[i_bin].mean = AVG(fg[a[ind_in_bin].master_index].color) -mnclr
        if ct ne 1 then color_list[i_bin].median = MEDIAN(fg[a[ind_in_bin].master_index].color)-mdclr else color_list[i_bin].median = fg[a[ind_in_bin].master_index].color-mdclr 
        color_list[i_bin].mean_err = STDDEV(a[ind_in_bin].color)/sqrt(n_elements(ind_in_bin))
        
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
	

    x = r_vector.mean_2d*1000
    y = color_list.mean
    y = y; - y[n_r_bin-1]
    y_err = color_list.mean_err
    print,y
    print,color_list.count
	th = 5
    my_Yr = [1e-5,1e-1]
    my_xr = [0.02,3.0]*1000.
	cc = ['g-r', 'g-W1', 'g-W2']
    if keyword_set(ps) then psopen, datapath+'reddening_'+cc[wise]+'m' + string(mmin, f='(F4.1)') +'--'+ string(mmax, f='(F4.1)') + 's' + string(smin*(-1), f='(F4.1)') +'--'+ string(smax*(-1), f='(F4.1)'), /helvetica, xsi=9, ysi=6, /inches, /color, /encapsulated
    loadct, 0
    if ps then !p.font=0 else !p.font = (-1)
    plot,x,y,/xlog,psym=4,yr=my_yr,ylog=1,xr=my_xr,$
      xtit=textoidl('separation [kpc]'),$
      ytit=my_y_tit, /xs, thick=th, xthick=th, ythick=th, /nodata
;    my_oploterr,x,y,y_err,psym=4,miny=my_yr[0], thick=th
    y = color_list.median
    print,y
    oplot,x*1.05,y,psym=4, color=getcolor('red',1)
        my_oploterr,x*1.05,y,y_err,psym=4,miny=my_yr[0],errcolor=getcolor('red',1)
    ylast = y[n_elements(y)-1]
  ;  oplot,x*1.05,y-ylast,psym=4, color=getcolor('green',1)
   ; 	my_oploterr,x*1.05,y-ylast,y_err,psym=4,miny=my_yr[0],errcolor=getcolor('green',1)
	
	xax = alog10(findgen(100)*10)
	oplot, 10^xax, 0.5*4.14d-3/3.1*((10.^xax)/100.)^(-0.84), color=200, thick=2*th, lines=1
	oplot, 10^xax, 4.14d-3/3.1*((10.^xax)/100.)^(-0.84), color=200, thick=2*th, lines=2
	oplot, 10^xax, 5*4.14d-3/3.1*((10.^xax)/100.)^(-0.84), color=200, thick=2*th, lines=1
	xyouts, 0.7, 0.8, mean(10^(a.mass_target)), /norm, charsize=3-ps*2
	xyouts, 0.7, 0.7, mean(a.z_target), /norm, charsize=3-ps*2
	print, 	(mean(10^(a.mass_target)))
	if keyword_set(ps) then psclose
endif


end
