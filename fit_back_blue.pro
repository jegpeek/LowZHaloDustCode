; fit the back bluening trend

dm = 0.1
minmag = 18
fts = fltarr(4, 31)
ftmags = fltarr(31)
for i=0, 30 do begin
	 measure_reddening, 1, ft, 0, /dofit, /gw, zbuf=0.015,/angle, /back, dmax= minmag+i*dm+2, dmin=minmag+i*dm, phmag=phmag, /dohist
	 ftmags[i] = phmag
	 fts[*, i] = ft
endfor


end