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

## Graphics modules
import matplotlib.pyplot as plt
import datetime
import numpy as np

## COMMAND-LINE ARG MODULES
import getopt
import argparse
import collections

## ASTRO
from astropy.io import fits
import scipy.stats as stats
#from scipy.stats import kurtosis, skew
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats 
##################################################



#### GET SCRIPT ARGS ####
def get_args():
	"""This function parses and return arguments passed in"""
	parser = argparse.ArgumentParser(description="Parse args.")
	parser.add_argument('-i', '--file', dest='filename', required=True, type=str, action='store',help='file name')
	
	args = parser.parse_args()	

	return args

##############
##   MAIN   ##
##############
def main():
	"""Main function"""

	## Get script args
	print('Get script args')
	try:
		args= get_args()
	except Exception as ex:
		print("Failed to get and parse options (err=%s)",str(ex))
		return 1

	filename= args.filename

	print("*** ARGS ***")
	print("filename: %s" % filename)
	print("************")
		

	## Read FITS	
	hdu= fits.open(filename, memmap=False)
	img = hdu[0].data

	## Get pixel list
	print("Getting pixel list ...")
	x= np.ravel(img)
	print(x)

	## Compute stats
	npixels= np.size(x)
	pixel_min= np.min(x)
	pixel_max= np.max(x)
	mean= np.mean(x)
	stddev= np.std(x,ddof=1)
	median= np.median(x)
	mad= mad_std(x)
	skewness= stats.skew(x)
	kurtosis= stats.kurtosis(x)
	print 'n=',npixels
	print 'min/max=',pixel_min,'/',pixel_max
	print 'mean=',mean
	print 'stddev=',stddev
	print 'median=',median
	print 'mad=',mad
	print 'skew=',skewness
	print 'kurtosis=',kurtosis

	## Compute robust stats
	niter= 1
	sigmaclip= 3
	[mean_clipped, median_clipped, stddev_clipped] = sigma_clipped_stats(x, sigma=sigmaclip, iters=niter, std_ddof=1)

	print 'mean_clipped=',mean_clipped
	print 'median_clipped=',median_clipped
	print 'stddev_clipped=',stddev_clipped

	
			
###################
##   MAIN EXEC   ##
###################
if __name__ == "__main__":
	#main()
	sys.exit(main())
