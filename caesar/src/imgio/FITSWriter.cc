// ***********************************************************************
// * License and Disclaimer                                              *
// *                                                                     *
// * Copyright 2016 Simone Riggi																			   *
// *																																	   *
// * This file is part of Caesar. 																		   *
// * Caesar is free software: you can redistribute it and/or modify it   *
// * under the terms of the GNU General Public License as published by   *
// * the Free Software Foundation, either * version 3 of the License,    *
// * or (at your option) any later version.                              *
// * Caesar is distributed in the hope that it will be useful, but 			 *
// * WITHOUT ANY WARRANTY; without even the implied warranty of          * 
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                *
// * See the GNU General Public License for more details. You should     * 
// * have received a copy of the GNU General Public License along with   * 
// * Caesar. If not, see http://www.gnu.org/licenses/.                   *
// ***********************************************************************
/**
* @file FITSWriter.cc
* @class FITSWriter
* @brief FITSWriter
*
* FITS Image Writer class
* @author S. Riggi
* @date 20/01/2015
*/

#include <FITSWriter.h>
#include <Img.h>

//Python interface
#include <TPython.h>
#include <Python.h>
//#include <numpy/ndarrayobject.h>
//#include <numpy/ndarraytypes.h>

#include <TObject.h>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <math.h>
using namespace std;

ClassImp(Caesar::FITSWriter)

namespace Caesar {


FITSWriter::FITSWriter() {

}//close costructor


FITSWriter::~FITSWriter() {

}//close destructor

int FITSWriter::Init(){
	
	try {
		//Initialize
		//Py_Initialize();

		//## Import python modules
		//== Numpy ==
		if(!TPython::Exec("import numpy as np")){
			cerr<<"FITSWriter::Init(): ERROR: Failed to import numpy python module!"<<endl;
			return -1;
		}
		//== pyfits ==
		if(!TPython::Exec("import pyfits")){
			cerr<<"FITSWriter::Init(): ERROR: Failed to import pyfits python module!"<<endl;
			return -1;
		}
	}
	catch(...){
		cerr<<"FITSWriter::Init(): ERROR: Failed to initialize python interface!"<<endl;
		return -1;
	}

	return 0;

}//close Init()


int FITSWriter::WriteFITS(Img* img,std::string outfilename){

	//## Check image
	if(!img){
		cerr<<"FITSWriter::WriteFITS(): ERROR: Null ptr to given image!"<<endl;
		return -1;
	}
	if(outfilename==""){
		cerr<<"FITSWriter::WriteFITS(): ERROR: Empty output filename!"<<endl;
		return -1;
	}

	//## Initialize and load needed python modules
	if(Init()<0){
		cerr<<"FITSWriter::WriteFITS(): ERROR: Initialization failed!"<<endl;
		return -1;	
	}

	//## Load data to python
	//outfilename
	PyObject* pystr = TPython::ObjectProxy_FromVoidPtr(&outfilename, "std::string");

	//image data
	int Nx= img->GetNbinsX();
	int Ny= img->GetNbinsY();
	std::vector< std::vector<float> > v;
	for(int j=0;j<Ny;j++){
		v.push_back( std::vector<float>() );
		for(int i=0;i<Nx;i++){
			double w= img->GetBinContent(i+1,j+1);
			v[j].push_back(w);
		}
	}
	PyObject* pyvec2d= TPython::ObjectProxy_FromVoidPtr(&v,"std::vector< std::vector<float> >");

	PyObject* pymain = PyImport_ImportModule("__main__");
	PyModule_AddObject(pymain, "outfile", pystr);	
	PyModule_AddObject(pymain, "mat", pyvec2d);
	Py_DECREF(pymain);
	
	//## Print the imported data
	if(!TPython::Exec("print 'outfile:',outfile")){
		cerr<<"FITSWriter::WriteFITS(): ERROR: Failed to import outfilename in python!"<<endl;
		return -1;
	}
	if(!TPython::Exec("print 'mat: ',mat")){
		cerr<<"FITSWriter::WriteFITS(): ERROR: Failed to import mat in python!"<<endl;
		return -1;
	}
	
  	
	//## Bind ROOT img obj
	cout<<"FITSWriter::WriteFITS(): INFO: Binding image to python..."<<endl;
	if(!TPython::Bind(img,"img")){
		cerr<<"FITSWriter::WriteFITS(): ERROR: Failed to bind image on the ROOT-Python interface!"<<endl;
		return -1;	
	}
	bool hasMetaData= img->HasMetaData();
	Caesar::ImgMetaData* metadata= 0;
	if(hasMetaData){
		cout<<"FITSWriter::WriteFITS(): INFO: Binding image meta-data to python..."<<endl;
		metadata= img->GetMetaData();
		if(!TPython::Bind(metadata,"metadata")){
			cerr<<"FITSWriter::WriteFITS(): ERROR: Failed to bind image metadata on the ROOT-Python interface!"<<endl;
			hasMetaData= false;
		}
	}

	
	//## Create a HDU object 	
	cout<<"FITSWriter::WriteFITS(): INFO: Creating a HDU from the input image..."<<endl;
	if(!TPython::Exec("hdu = pyfits.PrimaryHDU(mat)")){
		cerr<<"FITSWriter::WriteFITS(): ERROR: Failed to create hdu from bin array!"<<endl;
		return -1;
	}

	//## Fill header info	
	cout<<"FITSWriter::WriteFITS(): INFO: Filling hdu header..."<<endl;
	try {
		TPython::Exec("hdu.header[\'NAXIS1\']= metadata.Nx"); 
		TPython::Exec("hdu.header[\'NAXIS2\']= metadata.Ny"); 
		TPython::Exec("hdu.header[\'CTYPE1\']= metadata.CoordTypeX"); 
		TPython::Exec("hdu.header[\'CTYPE2\']= metadata.CoordTypeY"); 
		TPython::Exec("hdu.header[\'CRPIX1\']= metadata.Cx"); 
		TPython::Exec("hdu.header[\'CRPIX2\']= metadata.Cy"); 
		TPython::Exec("hdu.header[\'CRVAL1\']= metadata.Xc"); 
		TPython::Exec("hdu.header[\'CRVAL2\']= metadata.Yc"); 
		TPython::Exec("hdu.header[\'CDELT1\']= metadata.dX"); 
		TPython::Exec("hdu.header[\'CDELT2\']= metadata.dY"); 
		TPython::Exec("hdu.header[\'CROTA1\']= metadata.RotX"); 
		TPython::Exec("hdu.header[\'CROTA2\']= metadata.RotY"); 
		TPython::Exec("hdu.header[\'BUNIT\']= metadata.BUnit"); 
		TPython::Exec("hdu.header[\'BMAJ\']= metadata.Bmaj"); 
		TPython::Exec("hdu.header[\'BMIN\']= metadata.Bmin"); 
		TPython::Exec("hdu.header[\'BPA\']= metadata.Bpa"); 
		TPython::Exec("hdu.header[\'EPOCH\']= metadata.Epoch"); 
	}//close try
	catch(...){
		cerr<<"FITSWriter::WriteFITS(): ERROR: Failed to fill hdu header vars!"<<endl;
		return -1;
	}

	//## Write to FITS file
	cout<<"FITSWriter::WriteFITS(): INFO: Writing to fits..."<<endl;
	if(!TPython::Exec("hdu.writeto(str(outfile))")){
		cerr<<"FITSWriter::WriteFITS(): ERROR: Failed to write FITS file!"<<endl;
		return -1;
	}


	//## Finalize
	cout<<"FITSWriter::WriteFITS(): INFO: Finalizing..."<<endl;
	//Py_DECREF(pymain);
	//Py_Finalize();
	
	return 0;

}//close WriteFITS()




}//close namespace
