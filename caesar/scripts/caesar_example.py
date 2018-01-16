#!/usr/bin/python

##########################
## Import modules
#########################
import os, sys
import numpy as np
import pyfits
from ROOT import gSystem
#gSystem.Load('/home/riggi/Analysis/SKAProjects/SKATools/caesar-install/lib/libCaesar')
gSystem.Load('libCaesar')

from ROOT import Caesar
#######################

Nx= long(10)
Ny= long(20)
imgname= 'img'
xlow= float(0.)
ylow= float(0.)

try:
	img= Caesar.Image(Nx,Ny,xlow,ylow,imgname)
except Exception as e: 
	print(e)

Nx= img.GetNx()
Ny= img.GetNy()

print 'Nx=',Nx
