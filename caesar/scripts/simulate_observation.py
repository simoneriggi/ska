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

## ASTRO
#from scipy import ndimage
#import pyfits
#from astropy.io import fits
#from astropy.units import Quantity
#from astropy.modeling.parameters import Parameter
#from astropy.modeling.core import Fittable2DModel
#from astropy.modeling.models import Box2D, Gaussian2D, Ring2D, Ellipse2D, TrapezoidDisk2D, Disk2D, AiryDisk2D
#from photutils.datasets import make_noise_image

## ROOT
#import ROOT
#from ROOT import gSystem, TFile, TTree, gROOT, AddressOf

## CAESAR
#gSystem.Load('libCaesar')
#from ROOT import Caesar

## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections

## Graphics modules
import matplotlib.pyplot as plt
import pylab
##################################################


##############################
##   SIMULATE OBSERVATION   ##
##############################
def simulate_observation(skymodel_file,output_file,exec_simobs_step=True,skymodelmask_file='',total_time='43200s',telconfigs=['atca_6a.cfg','atca_6b.cfg','atca_ew352.cfg','atca_ew367.cfg'],imgaxes=['254.851041667','-41.4765888889','2.1GHz','I'],niter=100,interactive=False):
	"""Simulate an observation from a sky model"""

	#===========================
	#==   Print args
	#===========================
	print("*** ARGS ***")
	print("Input model filename: %s " % skymodel_file)
	print("Output filename: %s " % output_file)
	print("Total time (s): %s " % total_time)
	print("Tel configs: ", telconfigs)	
	print("************")

	## Strip filename and get basename
	skymodel_file_base, skymodel_file_ext= os.path.splitext(os.path.basename(skymodel_file))


	## Import sky model FITS file
	skymodel_img= str(skymodel_file_base) + '-casa'
	print ('Importing sky model FITS file %s in CASA as image %s...' % (skymodel_file,skymodel_img))
	importfits(fitsimage=skymodel_file,imagename=skymodel_img,defaultaxes=True,defaultaxesvalues=imgaxes,overwrite=True)

	## If given import sky model mask FITS image
	#skymodelmask_img= ''
	#if skymodelmask_file:
	#	skymodelmask_file_base, skymodelmask_file_ext= os.path.splitext(os.path.basename(skymodelmask_file))
	#	skymodelmask_img= str(skymodelmask_file_base) + '-casa'
	#	print ('Importing sky model mask FITS in CASA as %s ...' % skymodelmask_img)	
	#	importfits(skymodelmask_file,skymodelmask_img,defaultaxes=True,defaultaxesvalues=imgaxes,overwrite=True)


	## Setting simobserve parameters	
	project_name= 'simobs_' + skymodel_file_base

	## Generate simulated sky map visibilities
	print ('INFO: Generate simulated observation from sky model for given telescope configurations (n=%d configs)...' % len(telconfigs))
	vis_list= []	
	for config in telconfigs:
		config_base, config_ext= os.path.splitext(os.path.basename(config)) 
		vis= project_name + '.' + config_base + '.noisy.ms' 
		vis_list.append(vis) 
		if exec_simobs_step:
			print ('INFO: Simulating observation for tel config %s (vis=%s)' % (str(config),str(vis)))
			simobserve(project=project_name,complist='',skymodel=skymodel_img,direction='',antennalist=config,totaltime=total_time,obsmode='int',maptype='hexagonal',mapsize='',overwrite=True)

	vis_images= ''
	for i in range(0,len(vis_list)-1):
		vis_images+= vis_list[i] + ','
	vis_images+= vis_list[len(vis_list)-1]
	print ('INFO: visibilities=%s' % vis_images)

	## Cleaning
	print ('INFO: Cleaning simulated data and produce final image...')
	simanalyze(project=project_name,vis=vis_images, skymodel=skymodel_img,niter=niter, mask=skymodelmask_file, overwrite=True, interactive=interactive, analyze=True)	



