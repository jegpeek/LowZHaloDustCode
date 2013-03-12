function h2d_ri, xvec, yvec, dx, dy, zvec=zvec, xrng = xrng, yrng=yrng, zimg=zimg

	if ~keyword_set(xrng) then xrng= minmax(xvec)
	if ~keyword_set(yrng) then yrng= minmax(yvec)

	xvec = (xvec < (xrng[1]-1d-6)) > xrng[0]
	yvec = (yvec < (yrng[1]-1d-6)) > yrng[0]
	
	yvecscl = floor((yvec-yrng[0])/dy)
	xvecscl = floor((xvec - xrng[0])/dx)
	binsx = (xrng[1]-xrng[0])/dx
	binsy = (yrng[1]-yrng[0])/dy
	onedee = xvecscl + yvecscl*binsx
	
	hg = histogram(onedee, bin=1, nbins=binsx*binsy, min=0, reverse=ri)

	if keyword_set(zvec) then begin
		zimg = hg*0.0
		
		for i=0, binsx*binsy-1 do begin
			if (ri[i] ne ri[i+1]) then begin
				 zimg[i] = mean(zvec(ri[ri[i]:ri[i+1]-1]))
			endif
		endfor
		zimg = reform(zimg, binsx, binsy)
	endif
	hg = reform(hg, binsx, binsy)
	
	return, hg
end