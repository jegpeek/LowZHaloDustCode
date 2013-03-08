#!/usr/bin/env python

"""
"""
from termcolor import colored
import stomp,math,pyfits,sys,os,time,socket,numpy
import load_sample
import my_progress_bar
reload(load_sample)



############################################################
#### User defined parameters
############################################################

## Foreground (=target) file:
datapath = os.path.expanduser('~') + '/Dropbox/LowZHaloDustData/'
target_file = datapath + 'fg_MPAJHU.fits'

background_file = datapath + 'pg10.fits'
output_file_name = datapath + 'MPA-SDSS.fit'
output_file_name = datapath + 'MPA-SDSS_REVERSE.fit'

#background_file = datapath + 'g-W1_nod5.fits'
#output_file_name = datapath + 'MPA-WISE.fit'
#output_file_name = datapath + 'MPA-WISE_REVERSE.fit'


# a stomp-specific format. Created by Ryan. Code exists for making such masks within stomp, fwiw
## Map file (can be skipped)
#data_path = 'raid-cita/menard/DATA/SDSS/Galaxy/photo/DR7/photoz/eta.physics.ucdavis.edu/DR7/'
#map_file = data_path + 'stripe_photoz.hmap_basic'

#### Parameters
r_p_min_kpc = 20.
r_p_max_kpc = 3000.
theta_min_arcsec = 2.
delta_z = 0.0033
z_min = 0.04
############################################################
############################################################


print "Current output file: %s " %output_file_name
print "New output? [y/n]"
choice = sys.stdin.readline()
choice
if choice[0] == 'y':
	print "new output filename:"
	output_file_name = sys.stdin.readline().replace('\n','')
	print colored("Output filename used: %s",'cyan') %output_file_name
	
#    print 'Loading map: %s' %map_file
#    my_map = stomp.Map(map_file)
#    print colored("Map area: %3.2f deg^2" %my_map.Area() , 'cyan' )

print "Load foreground data? [1=Yes]"
choice = sys.stdin.readline()
if choice[0] == '1':
        target_sample, target_map = load_sample.load_sample(target_file,z_min=z_min) #,n_max=20000)
		
print "Load background data? [1=Yes]"
choice = sys.stdin.readline()
if choice[0] == '1':
	background_sample, background_map = load_sample.load_sample(background_file,z_min=z_min) #, n_max=20000)
	
print "Median foreground redshift:",numpy.median(target_sample.field('z'))
	
##############################
# A quad tree to store the positions of the background sources -- the things you look for.
## We create the tree
background_tree_map = stomp.IndexedTreeMap()
for background in background_map:
	background_tree_map.AddPoint(background)
	
	
##############################
# Initializing a bunch of parameters
indices = stomp.IndexVector()
neighbor_list = []
separation_list = []
physical_separation_list = []
color_list = []
z_target_list = []
mass_target_list = []
ssfr_target_list = []
z_background_list = []
master_index_list = []
master_list = []


print 'Finding pairs and computing physical separation...'
i_target = -1L
n_target = target_map.size()
my_progress = my_progress_bar.progressBar(0,n_target, 50)
# loop over each target (thing you sit on)
for target in target_map[:]:
	# target here is an element of the array of objects target_map
	my_progress.updateAmount_and_write_if_needed(i_target)
	i_target += 1
	
	z_target = target_sample[target.Index()].field('z')
	mass_target = target_sample[target.Index()].field('mass')
	ssfr_target = target_sample[target.Index()].field('ssfr')
	dist_to_target_Mpc = stomp.Cosmology_AngularDiameterDistance( numpy.double(z_target) )
	
	theta_min = r_p_min_kpc / (1e3 * dist_to_target_Mpc) * stomp.RadToDeg
	theta_max = r_p_max_kpc / (1e3 * dist_to_target_Mpc) * stomp.RadToDeg
	angular_bin_temp = stomp.AngularBin(theta_min,theta_max)
	
	# Find pairs is the brains of the code, and is a method of the tree object
	background_tree_map.FindPairs(target, angular_bin_temp, indices)


        # now let's look through all the objects...
	for index in indices:
		z_background = background_sample[index].field('z')
		
		# so to implement a redshift cut by target we do the tree search and then the cut. Makes sense.
		################### WARNING ###################
		# THIS IS THE BACKWARD HACK!!!!!!!!!!!!!!!
		# LOOKING FOR REDDENING BY GALAXIES BEHIND THE BACKGROUND
		# IS SIMPLY A CHECK TO MAKE SURE WE AREN'T BEING FOOLED
		################### WARNING ###################
		if ( z_background < (z_target - delta_z)):
			x = background_sample[index].field('RA')
			y = background_sample[index].field('DEC')
			tmp_ang = stomp.AngularCoordinate(numpy.double(x),numpy.double(y),stomp.AngularCoordinate.Equatorial)
			separation_arcsec = stomp.AngularCoordinate.AngularDistance(target, tmp_ang) * 3600.
			
			if (separation_arcsec > theta_min_arcsec):
				# append all the relevant information
				master_index_list.append(index)
				separation_list.append( separation_arcsec )
				z_target_list.append(z_target)
				mass_target_list.append(mass_target)
				ssfr_target_list.append(ssfr_target)
				z_background_list.append(z_background)
				physical_separation = stomp.Cosmology_ProjectedDistance( numpy.double(z_target), separation_arcsec/3600.)
				# the physical separation is in Mpc
				physical_separation_list.append(physical_separation)
				color = background_sample[index].field('color')
				color_list.append(color)
				
				master_list.append([z_target,mass_target, ssfr_target, index,separation_arcsec,physical_separation,color])

				
print "found %i pairs" % len(z_background_list)
			

#write the fits file
index_col = pyfits.Column(name="master_index", format="J", array=master_index_list)
z_target_col = pyfits.Column(name="z_target", format="E", array=z_target_list)
mass_target_col = pyfits.Column(name="mass_target", format="E", array=mass_target_list)
ssfr_target_col = pyfits.Column(name="ssfr_target", format="E", array=ssfr_target_list)
z_background_col = pyfits.Column(name="z_background", format="E", array=z_background_list)
physical_separation_col = pyfits.Column(name="physical_separation_Mpc", format="E", array=physical_separation_list)
color_col = pyfits.Column(name="color", format="E", array=color_list)

cols = pyfits.ColDefs([index_col, z_target_col, mass_target_col, ssfr_target_col, z_background_col, physical_separation_col,color_col])

my_output = pyfits.new_table(cols)
print "Writing to %s..." % output_file_name
if os.path.exists(output_file_name):
	print 'Creating new file'
	os.remove(output_file_name)
my_output.writeto(output_file_name)

print "\a\a\a"
