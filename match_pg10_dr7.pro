; match pg10 dr7

datapath = '~/Dropbox/LowZHaloDustData/'

file2 = datapath + 'pg10.fits'
pg10 = MRDFITS(file2,1, /sil)

dr7 = mrdfits('~/Documents/halodust/gal_info_dr7_v5_2.fits', 1)
close_match_radec,pg10.ra,pg10.dec,dr7.ra,dr7.dec,m1,m2,2./3600.,1,miss1

save, m2, f='pg10_dr7_match.sav'

file_galaxy = datapath +'fg_MPAJHU.fits'



end