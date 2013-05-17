pro plottruncfit, p, mass_lstar, ssfr_11, opl=opl, _REF_EXTRA=_extra, y=y, xax=xax

p0 = reform(p)
xax = 10^((findgen(100)+1)/100*2.5 + 1); kpc
if n_elements(p0) eq 4 then p0 = [p0, 1e6]
y = p0[0]*(xax/100)^(p0[1])*mass_lstar^(p0[2])*ssfr_11^(p0[3])*exp((-1)*(xax/100)/p0[4])

if ~keyword_set(opl) then begin
	plot, xax, y, /xlog, /ylog, yra=[1d-4, 1d-1], /xs, xtitle='radius kpc', ytitle='E(B-V) equiv'
endif else begin 
	oplot, xax, y, _extra=_extra
endelse
end


