pro plot_utrends, R_fid, ps=ps

datapath = '~/Dropbox/LowZHaloDustData/'
vals = fltarr(5)
errs = fltarr(5, 2)
clr = [2, 4, 8, 9]
for i=0, 4 do begin
	if i eq 0 then restore, datapath + 'MC_WISE3_fwd.sav' else restore, datapath + 'MC_SDSSfwd'+string(clr[i-1], f='(I2.2)')+'.sav'
	compute_scatter, allp, R_fid, 1, 1, val, err
	vals[i] = val
	errs[i, *] = err
endfor

circle, /fill
plot, findgen(5), vals, syms=3, psym=8, xra=[-1, 5], /xs, yra=minmax(errs)
for i=0, 4 do oplot, [i, i], errs[i, *], thick=2

c1 = [1,2,0,3,0,1,1,2, 0, 0]
c2 = [2,3,1,4,2,3,4,4, 3, 4]
fc = [5.15500   ,   3.79300  ,    2.75100   ,   2.08600   ,   1.47900]

factors = [fc[0]-0.015, fc[0] -fc[1], fc[0] -fc[2], fc[0] -fc[3], fc[0] -fc[4]]

exts= fltarr(5)
exterrs = fltarr(5, 2)

exts[0] = vals[0]*factors[0]
exterrs[0, *] = errs[0, *]*factors[0]	
for i=1, 4 do exts[i] = exts[0] - vals[i]*factors[i]

ee0 = (abs(errs[0, 0]-vals[0]) +  abs(errs[0, 1]-vals[0]))/2.*factors[0]

for i=1, 4 do begin
	errfact= (abs(errs[i, 0]-vals[i]) +  abs(errs[i, 1]-vals[i]))/2.*factors[i]
	exterrs[i, *] = exts[i] + [-1, 1]*errfact;(sqrt(errfact^2 + ee0^2))
endfor
if keyword_set(ps) then psopen, 'red_law', /color, /landscape, /helvetica
!p.font = 0
th=5
plot, findgen(5), exts, syms=3, psym=8, xra=[-1, 5], /xs, yra=minmax(exterrs), /nodata, xtickn = [' ', 'u', 'g', 'r', 'i', 'z', ' '], ytitle='Extinction at '+string(R_fid, f='(I2.2)')+' kpc', xthick=th, ythick=th, thick=th

polyfill, [findgen(5), reverse(findgen(5))], [exts - exts*(exts[0]- exterrs[0,0])/exts[0], reverse(exts + exts*(exterrs[0,1]-exts[0])/exts[0])], color=200

oplot, findgen(5), exts, syms=3, psym=8, thick=th

for i=1, 4 do oplot, [i, i], exterrs[i, *], thick=th

;oplot, findgen(5), exts - exts*(exts[0]- exterrs[0,0])/exts[0]
;oplot, findgen(5), exts + exts*(exterrs[0,1]-exts[0])/exts[0]

oplot, findgen(5), fc/180., line=1, thick=th

if keyword_set(ps) then psclose

end