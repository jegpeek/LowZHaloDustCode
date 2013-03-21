; plot the number of galaxies within some physical impact parameter
readdata=0
if readdata eq 1 then begin
	datapath = '~/Dropbox/LowZHaloDustData/'
	
	afwd = mrdfits(datapath + 'MPA-SDSS_angspace.fit', 1)
	abck = mrdfits(datapath + 'MPA-SDSS_REVERSE_angspace.fit', 1)
	
;	restore, datapath +'MPA-SDSS_asec.sav'
;	as_f = asec
	
;	restore, datapath +'MPA-SDSS_REVERSE_asec.sav'
;	as_b = asec

	
	zt = [afwd.z_target, abck.z_target]
	zb = [afwd.z_background, abck.z_background]
endif

psopen, 'histo_zbuf', xsi=8.5, ysi=3.3, /inches, /helvetica, /encap
!p.font=0
th=2
plothist, zb-zt, bin=0.001, xra=[-0.25, 0.25], thick=th, xthick=th, ythick=th, /ylog, yra=[100, 1e6]
plothist, afwd.z_background-afwd.z_target, bin=0.001, /over, /fill, fcolor=170, /ylog
plothist, abck.z_background-abck.z_target, bin=0.001, /over, /fill, fcolor=90, /ylog

whoutfwd = where((afwd.z_background-afwd.z_target) lt 0.015)
whoutbck = where((abck.z_background-abck.z_target) gt (-0.015))

plothist, (afwd.z_background-afwd.z_target)[whoutfwd], bin=0.001, /over, /fill, fcolor=220, /ylog
plothist, (abck.z_background-abck.z_target)[whoutbck], bin=0.001, /over, /fill, fcolor=220, /ylog

psclose

end