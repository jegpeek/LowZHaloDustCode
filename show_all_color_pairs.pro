; show all pairs
psopen, 'reddening_10-color', /helvetica, /color, /encap
!p.font=0
!p.multi=[0, 2, 5]

for j=0, 9 do begin
 	measure_reddening, 0, ft, 0, /dohist, zbuf= 0.015, use10=0.1+j
endfor
psclose

psopen, 'reverse_10-color', /helvetica, /color, /encap
!p.font=0
!p.multi=[0, 2, 5]

for j=0, 9 do begin
 	measure_reddening, 0, ft, 0, /dohist, zbuf= 0.015, use10=0.1+j, /back
endfor
psclose



end