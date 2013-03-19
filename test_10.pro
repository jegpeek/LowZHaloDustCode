; test dz effects

getdata = 1
if getdata then begin
		
	ft10_bck = fltarr(4, 10, 5)
	ft10_fwd = fltarr(4, 10, 5)
	for j=0, 4 do begin

	for i=0, 9 do begin
		measure_reddening, 0, ft, 0, /dofit, a0=a0_10, zbuf=0.01*j, use10=i
		ft10_fwd[*, i, j] = ft 
		print, '########' + string(i)
	endfor

	for i=0, 9 do begin
		measure_reddening, 0, ft, 0, /dofit, /back, a0=a0_10bck, zbuf=0.01*j, use10=i
		ft10_bck[*, i, j] = ft 
		print, '########' + string(i)
	endfor
	endfor

endif



end