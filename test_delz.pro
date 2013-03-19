; test dz effects

getdata = 0
if getdata then begin
	
	for i=0, 49 do begin
		measure_reddening, 1, ft, 1, /dofit, a0=a0wise, zbuf=0.001*i 
		ft50_fwd[*, i] = ft 
		print, '########' + string(i)
	endfor

	for i=0, 49 do begin
		measure_reddening, 1, ft, 1, /dofit, /back, a0=a0wiseback, zbuf=0.001*i 
		ft50_bck[*, i] = ft 
		print, '########' + string(i)
	endfor

endif

zs = findgen(50)*0.001
clr = 100
!p.multi=[0, 1, 5]
loadct, 0
th=4
plothist, a2.z_target-a2.z_background, bin=0.001, xra=[-0.25, 0], charsize=2, thick=th
plothist, a5.z_target-a5.z_background, bin=0.001, xra=[0, 0.25], charsize=2, color=clr, thick=th


for j=0, 2 do begin
	plot, zs , ft50_fwd[j, *],  xtitle='zbuffer', thick=th, charsize=2
	oplot, zs, ft50_bck[j, *], thick=th,color=clr
endfor

end