; make new PG10 colors with fewer assumptions
datapath = '~/Dropbox/LowZHaloDustData/'

; some useful metadata
hvc = mrdfits(datapath + 'hvc_data_dr7_1sig.fits', 1, hdr)
dev = mrdfits(datapath + 'dev_dr7_1sig.fits', 1, hdr) 
img = mrdfits(datapath + 'imaging_dr7_1sig.fits', 1, hdr)  
 
c1 = [1,2,0,3,0,1,1,2, 0, 0]
c2 = [2,3,1,4,2,3,4,4, 3, 4]
clr = ['u', 'g', 'r', 'i', 'z']


; note we are assuming the standard SFD extinction values here:
 fc = [5.15500   ,   3.79300  ,    2.75100   ,   2.08600   ,   1.47900]

medsm = 500
ord = 4
ps = 0
th=3
!p.multi=[0, 10, 3, 1, 1]
cs = 2
model = 1

ngcs = fltarr(n_elements(hvc), 10)
rats = fltarr(10)
ogwid = [0.0243, 0.0307, 0.0708, 0.0384, 0.0476, 0.0253, 0.0247, 0.0276, 0.0424, 0.0379]

for i=0, 9 do begin

	if model eq 1 then begin
		wh = findgen(n_elements(img))
		nct = n_elements(img)
		ctwh = nct
		mmag1 = 22.5 - 2.5*alog10(img.modelflux[c1[i]])
		mmag2 = 22.5 - 2.5*alog10(img.modelflux[c2[i]])
		colori = (((mmag1-mmag2) -  (hvc.extinction[c1[i]]- hvc.extinction[c2[i]]))[wh])
	endif else begin
		wh = where(dev.devmag[c1[i]]-dev.devmag[c2[i]] ne 0, ctwh)
		nct = n_elements(wh)
		colori = (((dev.devmag[c1[i]]-dev.devmag[c2[i]]) -  (hvc.extinction[c1[i]]- hvc.extinction[c2[i]]))[wh])
	endelse
	nct = n_elements(wh)
	mwh = (findgen(nct))[medsm/2:nct-medsm/2-1]

	zmed = median(   colori[sort(hvc[wh].z)], medsm)
	loadct, 0, /sil 
	plot,  colori, hvc[wh].z, psym=3, ytitle='redshift' ,xra=[0, 5], /xs, charsize=cs, xtitle=clr[c1[i]] + '-'+ clr[c2[i]]
	loadct, 13, /sil
	oplot, zmed[mwh], ((hvc[wh].z)[sort(hvc[wh].z)])[mwh], color=ps*255, thick=th
	pf = poly_fit(((hvc[wh].z)[sort(hvc[wh].z)])[mwh], zmed[mwh], ord, yfit=yf)
 
	oplot, yf, (hvc[wh].z)[sort(hvc[wh].z)], color=fitcolor, thick=th
	residc = colori - total(rebin(reform(hvc[wh].z, ctwh, 1), ctwh, ord+1)^(rebin(reform(findgen(ord+1), 1, ord+1), ctwh, ord+1))*rebin(pf, ctwh, ord+1), 2)
	loadct, 0, /sil
	plot,  residc, hvc[wh].absmag[2], psym=3, ytitle='r-band absmag', xra=[-1, 1], /xs	, charsize=cs, xtitle='corrected' + clr[c1[i]] + '-'+ clr[c2[i]]
	plothist, residc < 2 > (-2), xh, yh, bin=0.001, xra=[-1, 1], charsize=cs
	fitg = gaussfit(xh, yh, gf, nterms=4)
	loadct, 13, /sil
	oplot, xh, fitg, thick=th
	loadct, 0, /sil
	xyouts, !x.crange[1]*0.35, !y.crange[1]*2./3., clr[c1[i]] + '-'+ clr[c2[i]], charsize=2
	xyouts, !x.crange[1]*0.35, !y.crange[1]/2., string(gf[2]*1000/(fc[c1[i]] - fc[c2[i]]), f='(I3.3)'), charsize=2
	ngcs[*, i] = residc/(fc[c1[i]] - fc[c2[i]])
	print, clr[c1[i]] + '-'+ clr[c2[i]]+ ': ' + string((gf[2]/(fc[c1[i]] - fc[c2[i]]))/ogwid[i])
	rats[i] = (gf[2]/(fc[c1[i]] - fc[c2[i]]))/ogwid[i]
endfor
 
save, ngcs, f='ngcs.sav'


end