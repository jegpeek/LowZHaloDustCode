PRO DISPLAY_DATA,galaxy,gW1,pg10

!p.charsize=1.5

!p.multi=[0,2,2]

n_gal_plot = 1e3
iseed = 1.2
ind_random = RANDOMU(iseed,n_gal_plot)*n_elements(galaxy)

plothist,galaxy[ind_random].z,bin=0.01,xr=[0,0.2]
plothist,gW1.z,bin=0.01,xtit='redshift'

my_xr = 0.4*[-1,1]
plothist,gW1.color,bin=0.01,xtit='g-W1',xr=my_xr
plothist,pg10.color,bin=0.01,xtit='(g-r)pg10',xr=my_xr

print,STDDEV(gW1.color)
print,STDDEV(pg10.color)
pause

!p.multi=0

n_gal_plot = 300
ind_random = RANDOMU(iseed,n_gal_plot)*n_elements(galaxy)

x = galaxy[ind_random].ra
y = galaxy[ind_random].dec
plot,x,y,psym=1

ind_random = RANDOMU(iseed,n_gal_plot)*n_elements(gW1)
x = gW1[ind_random].ra
y = gW1[ind_random].dec
oplot,x,y,psym=1,color=getcolor('red',1)

ind_random = RANDOMU(iseed,n_gal_plot)*n_elements(pg10)
x = pg10[ind_random].ra
y = pg10[ind_random].dec
oplot,x,y,psym=1,color=getcolor('blue',1)


end
