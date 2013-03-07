pro rollmed, x, y, dx, xcent, meds, nerr, dn=dn, xmm=xmm


if keyword_set(dn) then begin
	srty = y(sort(x))
	srtx = x(sort(x))
	;med = median(srty, dn)
	;wh = where(med ne srty)
	;xcent = (x(sort(x)))[wh]
	;meds = med[wh]
	nbin = floor(n_elements(y)/float(dn))
	meds = fltarr(nbin)
	xcent = fltarr(nbin)
	for i=0l, nbin-1 do begin
		meds[i] = median(srty[i*dn:(i+1)*dn-1])
		xcent[i] = median(srtx[i*dn:(i+1)*dn-1])
	endfor
endif else begin
if not keyword_set(xmm) then xmm = minmax(x)
nx = ceil((xmm[1]-xmm[0])/dx)
xcent = (findgen(nx)-0.5)*dx+xmm[0]
meds = fltarr(nx)
nerr = meds
for i=0l, nx-1l do begin
	wh = where((x ge xcent[i]-dx/2) and (x lt xcent[i]+dx/2), ct)
	if ct gt 10 then begin
		meds[i] = median(y[wh])
		plothist, y[wh], xw, yw, bin=0.1, /no
	;	rl = gaussfit(xw, yw, al, n=4)
	;	nerr[i] = al[2]/sqrt(float(ct))*sqrt(!pi/2.)
	endif
endfor
endelse
end
