#!/usr/bin/python

##################################################
###          MODULE IMPORT
##################################################
## STANDARD MODULES
import os
import sys
import subprocess
import string
import time
import signal
from threading import Thread
import datetime
import numpy as np
import random
import math

## ASTRO MODULES
from astropy.io import fits

def simulate_observation(concat_vis='merged_vis.ms',skymodel='skymodel.fits',exec_simobs_step=True,project_name='sim',total_time='43200s',telconfigs=['atca_6a.cfg','atca_6b.cfg','atca_ew352.cfg','atca_ew367.cfg'],obsmode='int',maptype='hexagonal',mapsize='',direction='',indirection='',incell='',incenter= '2.1GHz',inwidth='1GHz',integration='10s',import_model=True,imgaxes=['254.851041667','-41.4765888889','2.1GHz','I'],use_noise_vis=True,add_thermal_noise=False):
	"""Simulate an observation from a sky model"""

	#===========================
	#==   Print args
	#===========================
	print("*** ARGS ***")
	print("Input model filename: %s " % skymodel)
	print("Total time (s): %s " % total_time)
	print("Tel configs: ", telconfigs)	
	print("************")

	## Strip filename and get basename
	skymodel_file_base, skymodel_file_ext= os.path.splitext(os.path.basename(skymodel))

	## Import sky model FITS file?
	skymodel_img= str(skymodel_file_base) + '-casa'
	if import_model:	
		print ('Importing sky model FITS file %s in CASA as image %s...' % (skymodel,skymodel_img))
		###importfits(fitsimage=skymodel_file,imagename=skymodel_img,defaultaxes=True,defaultaxesvalues=imgaxes,overwrite=True)
		importfits(fitsimage=skymodel,imagename=skymodel_img,overwrite=True)

	## Generate simulated sky map visibilities
	print ('INFO: Generate simulated observation from sky model for given telescope configurations (n=%d configs)...' % len(telconfigs))
	vis_list= []	
	for config in telconfigs:
		config_base, config_ext= os.path.splitext(os.path.basename(config)) 
		if use_noise_vis:
			vis= project_name + '/' + project_name + '.' + config_base + '.noisy.ms'
		else:
			vis= project_name + '/' + project_name + '.' + config_base + '.ms'
		vis_list.append(vis) 

		# Execute the simulation step?
		if exec_simobs_step:

			# Add thermal noise?
			if add_thermal_noise:
				thermalnoise= 'tsys-atm'
			else:
				thermalnoise= ''

			# Simulate observation
			print ('INFO: Simulating observation for tel config %s (vis=%s)' % (str(config),str(vis)))
			#simobserve(project=project_name,skymodel=skymodel_img,direction='',antennalist=config,totaltime=total_time,obsmode='int',maptype='hexagonal',mapsize='',overwrite=True)
			simobserve(
				project=project_name, 
				skymodel=skymodel_img, 
				incenter=incenter, 
				incell=incell,
				inwidth=inwidth,
				indirection=indirection,
				direction=direction,	
				antennalist=config,
				totaltime=total_time,
				integration=integration,
				obsmode=obsmode,
				maptype=maptype,
				mapsize=mapsize,
				thermalnoise=thermalnoise,
				overwrite=True
			)

	## Concatenate visibilities
	print 'Concatenating visibilities: ', vis_list
	concat(vis=vis_list, concatvis=concat_vis)	





def analyze_observation(vis,project_name='rec',niter=500,mask='',skymodel='',fitsout='output.fits',interactive=False,analyze=True):
	""" Clean & anayze observation """
	## Cleaning
	print ('INFO: Cleaning simulated data and produce final image...')
	simanalyze(
		project=project_name,
		vis=vis, 
		skymodel=skymodel,
		niter=niter, 
		mask=mask,
		imdirection='', 
		interactive=interactive, 
		analyze=analyze,
		overwrite=True
	)	

	## Export to FITS
	vis_base, vis_ext= os.path.splitext(os.path.basename(vis)) 
	exported_map= project_name + '/' + project_name + vis_base + '.image' 
	print ('INFO: Exporting CASA map %s to FITS...' % exported_map)
	exportfits(imagename=exported_map, fitsimage=fitsout, history=False, overwrite=True)

	## Get image head & data
	#imghead= imhead(imagename=exported_map)
	#nx= imghead['shape'][0] # RA (x-axis)
	#ny= imghead['shape'][1] # DEC (y-axis)
	#boxselection= '0,0,' + str(nx) + ',' + str(ny)
	#imgdata= imval(imagename=exported_map,box=boxselection)

	## Set FITS data & header
	## NB: CASA exportfits uses the SIN projection in output, I want the same image format of the model! 
	##     Found that I need to transpose casa image data and reverse x-axis (fliplr)
	#data= np.fliplr(np.transpose(imgdata['data']))
	
	# Define FITS header
	#header= fits.Header()	
	#header.set('SIMPLE','T')
	#header.set('BITPIX','-32')
	#header.set('NAXIS1', str(self.nx))
	#header.set('NAXIS2', str(self.ny))
	#header.set('NAXIS3', 1)
	#header.set('NAXIS4', 1)
	#header.set('BUNIT', str(imghead['unit']))
	#header.set('BMAJ', imghead['restoringbeam']['major']['value'])
	#header.set('BMIN', imghead['restoringbeam']['major']['value'])
	#header.set('BPA', imghead['restoringbeam']['positionangle']['value'])
	#header.set('BSCALE',1.)
	#header.set('BZERO',0.)
	#header.set('CDELT1',self.pixsize/3600.)
	#header.set('CDELT2',self.pixsize/3600.)
	#header.set('CTYPE1',self.ctype1)
	#header.set('CTYPE2',self.ctype2)
	#header.set('CRPIX1',self.crpix1)
	#header.set('CRPIX2',self.crpix2)
	#header.set('CRVAL1',self.crval1)
	#header.set('CRVAL2',self.crval2)

	## Exporting to FITS
	#hdu = fits.PrimaryHDU(data=data,header=header)
	#hdulist = fits.HDUList([hdu])
	#hdulist.writeto(fitsout,overwrite=True)


def clean_observation(vis,niter=500,mask='',imagename='',field='',image_size='',cell_size='1arcsec',phase_center='',weighting='briggs',projection='CAR',deconvolver='clark',gridder='mosaic',fitsout='output.fits',interactive=False):
	""" Clean observation """
	
	## Set image name if empty
	if not imagename:
		project_name= 'rec'
		recimg= 'recmap'
		imagename=project_name + '/' + recimg	

	## Cleaning
	print ('INFO: Cleaning simulated data and produce final image...')
	tclean(
		imagename=imagename,
		vis=vis,
		field=field, 
		mask=mask,
		imsize=image_size,
		phasecenter=phase_center,
		projection=projection,
		cell=cell_size,	
		niter=niter,
		deconvolver=deconvolver,
		gridder=gridder,
		weighting=weighting,
		interactive=interactive
	)

	## Exporting to FITS
	exported_map= imagename + '.image' 
	print ('INFO: Exporting CASA map %s to FITS...' % exported_map)
	exportfits(imagename=exported_map, fitsimage=fitsout, history=False, overwrite=True)


def run_sim(outvis='vis.ms',skymodel='skymodel.fits',sim_projname='sim',total_time='43200s',telconfigs=['atca_all.cfg'],obsmode='int',maptype='hexagonal',mapsize='',direction='',indirection='',incell='',incenter= '2.1GHz',inwidth='1GHz',integration='10s',use_noise_vis=True,add_thermal_noise=False,niter=1000,recimgname='',mask=''):
	""" Simulate and imaging a sky model field"""

	## Simulate observation
	simulate_observation(
		concat_vis=outvis,
		skymodel=skymodel,
		project_name=sim_projname,
		total_time=total_time,
		telconfigs=telconfigs,
		obsmode=obsmode,
		maptype=maptype,
		mapsize=mapsize,
		direction=direction,
		indirection=indirection,
		incell=incell,
		incenter=incenter,
		inwidth=inwidth,
		integration=integration,
		use_noise_vis=use_noise_vis,
		add_thermal_noise=add_thermal_noise
	)

	# Retrieve number of pointings
	tb.open(outvis + '/FIELD')
	npointings= tb.nrows
	print ('#' + str(npointings) + ' pointings present'
	
	# Imaging
	clean_observation(
		vis=outvis,
		niter=niter,
		mask=mask,
		imagename=recimgname,
		#field='',image_size='',cell_size='1arcsec',phase_center='',weighting='briggs',projection='CAR',deconvolver='clark',gridder='mosaic',fitsout='output.fits',interactive=False
	)
	
	#lm=casac.linearmosaic()
	#lm.defineoutputimage(nx=2560, cellx='1arcsec', imagecenter='J2000 254.85067091580214deg -41.47631111052697deg', outputimage='linmos')
	#lm.setlinmostype('optimal')
	#lm.makemosaic(images=['rec_field0/recmap.image','rec_field1/recmap.image'],weightimages=['rec_field0/recmap.pb','rec_field1/recmap.pb'])
	#lm.saultweightimage(’test_sault.linmos’) 

def casa_to_fits(image_name,outfits='exported_map.fits'):
	""" Export a CASA image to FITS """

	## Get image head & data
	imghead= imhead(imagename=image_name,mode='list')
	nx= imghead['shape'][0] # RA (x-axis)
	ny= imghead['shape'][1] # DEC (y-axis)
	boxselection= '0,0,' + str(nx) + ',' + str(ny)
	imgdata= imval(imagename=image_name,box=boxselection)

	## Set FITS data & header
	## NB: CASA exportfits uses the SIN projection in output, I want the same image format of the model! 
	##     Found that I need to transpose casa image data and reverse x-axis (fliplr)
	data= np.fliplr(np.transpose(imgdata['data']))
	
	# Define FITS header
	header= fits.Header()	
	header.set('SIMPLE','T')
	header.set('BITPIX','-32')
	header.set('NAXIS1', str(nx))
	header.set('NAXIS2', str(ny))
	#header.set('NAXIS3', 1)
	#header.set('NAXIS4', 1)
	header.set('BUNIT', str(imghead['unit']))
	header.set('BMAJ', imghead['restoringbeam']['major']['value'])
	header.set('BMIN', imghead['restoringbeam']['major']['value'])
	header.set('BPA', imghead['restoringbeam']['positionangle']['value'])
	header.set('BSCALE',1.)
	header.set('BZERO',0.)
	#header.set('CDELT1',self.pixsize/3600.)
	#header.set('CDELT2',self.pixsize/3600.)
	#header.set('CTYPE1',self.ctype1)
	#header.set('CTYPE2',self.ctype2)
	#header.set('CRPIX1',self.crpix1)
	#header.set('CRPIX2',self.crpix2)
	#header.set('CRVAL1',self.crval1)
	#header.set('CRVAL2',self.crval2)

	## Exporting to FITS
	#hdu = fits.PrimaryHDU(data=data,header=header)
	#hdulist = fits.HDUList([hdu])
	#hdulist.writeto(fitsout,overwrite=True)

def simulate_observation_light():
	""" Simulate observation """	
	
	## Define options
	project_name= 'simobs'
	skymodel= 'skymodel.fits'
	skymodel_img= 'skymodel'
	mask= 'mask.txt'
	output_map= 'output.fits'
	incenter= '1.3GHz'
	inwidth= '1GHz'
	phase_center= 'J2000 254.851041667deg -41.4765888889deg'
	cell_size='1arcsec'
	projection='NCP'
	image_size= 1000
	##indirection='J2000 16h59m24 -41d28m34'
	indirection=''
	direction= ''
	tel_config= 'atca_6a.cfg'
	total_time='43200s'
	obsmode='int'
	maptype='square'
	mapsize= '' ## Same model size
	#imgaxes=['254.851041667','-41.4765888889','2.1GHz','I']
	imgaxes= ''
	niter= 500

	## Import sky model FITS file
	print ('INFO: Importing sky model FITS file %s in CASA as image %s...' % (skymodel,skymodel_img))
	#importfits(skymodel,skymodel_img,defaultaxes=True,defaultaxesvalues=imgaxes,overwrite=True)
	importfits(skymodel,skymodel_img,overwrite=True)

	## Simulate observation
	print ('INFO: Generate simulated observation from sky model for given telescope configurations')
	simobserve(
		project=project_name, 
		skymodel=skymodel_img, 
		incenter=incenter, 
		incell='',
		inwidth=inwidth,
		indirection=indirection,
		direction=direction,	
		antennalist=tel_config,
		totaltime=total_time,
		obsmode=obsmode,
		maptype=maptype,
		mapsize=mapsize,
		overwrite=True
	)

	## Cleaning
	#print ('INFO: Cleaning simulated data and produce final image...')
	#config_base, config_ext= os.path.splitext(os.path.basename(tel_config)) 
	#vis= project_name + '/' + project_name + '.' + config_base + '.noisy.ms' 
	#simanalyze(
	#	project=project_name,
	#	vis=vis, 
	#	skymodel=skymodel_img,
	#	niter=niter, 
	#	mask=mask,
	#	imdirection='', 
	#	overwrite=True, 
	#	interactive=False, 
	#	analyze=True	
	#)

	#cleaned_map= project_name + '/' + project_name + '.' + config_base + '.noisy.image' 
	#tclean(
	#	imagename=cleaned_map,
	#	vis=vis, 
	#	#mask=mask,
	#	imsize=image_size,
	#	#phasecenter=phase_center,
	#	projection=projection,
	#	cell=cell_size,	
	#	niter=niter,
	#	#deconvolver='clark',
	#	interactive=False
	#)

	## Write fits file
	#exported_map= project_name + '/' + project_name + '.' + config_base + '.noisy.image' 
	#exportfits(imagename=cleaned_map, fitsimage=output_map, history=False, overwrite=True)


