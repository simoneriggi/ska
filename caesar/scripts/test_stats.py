#!/usr/bin/python

##################################################
###          MODULE IMPORT
##################################################
import numpy as np
import scipy.stats as stats
##import astropy.stats as astrostats
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats 

## Generate data 
np.random.seed(0)
x = np.arange(200)
y = np.zeros(200)
c = stats.bernoulli.rvs(0.35, size=x.shape)
y += (np.random.normal(0., 0.2, x.shape) + c*np.random.normal(3.0, 5.0, x.shape))

print(y)

## Compute stats
n= np.size(x)
mean= np.mean(y)
stddev= np.std(y,ddof=1)
median= np.median(y)
mad= mad_std(y)
M1= stats.moment(y,1)
M2= stats.moment(y,2)*n
M3= stats.moment(y,3)*n
M4= stats.moment(y,4)*n
print 'mean=',mean
print 'stddev=',stddev
print 'median=',median
print 'mad=',mad
print 'M1=',M1
print 'M2=',M2
print 'M3=',M3
print 'M4=',M4


## Compute robust stats
niter= 1
sigmaclip= 3
[mean_clipped, median_clipped, stddev_clipped] = sigma_clipped_stats(y, sigma=sigmaclip, iters=niter, std_ddof=1)

print 'mean_clipped=',mean_clipped
print 'median_clipped=',median_clipped
print 'stddev_clipped=',stddev_clipped

