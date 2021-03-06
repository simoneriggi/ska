# CAESAR

Compact And Extended Source Automated Recognition

## **About**  
CAESAR is a C++ software tool for automated source finding in astronomical maps. 
It is distributed under the GNU General Public License v3.0. 
If you use CAESAR for your research, please acknowledge it in your papers by citing the following paper:

* S. Riggi et al., "Automated detection of extended sources in radio maps:
progress from the SCORPIO survey", MNRAS (2016) doi: 10.1093/mnras/stw982, arXiv:1605.01852

## **Status**
This is the CAESAR version used in the reference paper.  
**NB: Please note this is an old version. A new CAESAR version is under development (see https://github.com/simoneriggi/ska.git, caesar directory)**. Some of the functionalities, like extended source search, are yet to be integrated in the new version. The planned major improvements for the new version:

* Build, test and integration system
* Code refactoring and organization
* Logging system
* Parallel and optimized versions of extraction algorithms
* Removed dependency from VLFEAT library
* Python interface
* TBD 

## **Installation**  
Install the project dependencies:  
* ROOT [https://root.cern.ch/], to be built with FITSIO, PyROOT, RInterface options enabled. Make sure that the FindROOT.cmake is present in $ROOTSYS/etc/cmake directory after installation.
* OpenCV [http://opencv.org/]
* R [https://www.r-project.org/], install also these additional packages: RInside, Rcpp, rrcovHD, truncnorm, FNN, akima
* vlfeat [http://www.vlfeat.org/]
* boost [http://www.boost.org/]
* python (>=2.7) [https://www.python.org/], install also these additional modules: pyfits
* cmake (>=2.8) [https://cmake.org]  
  
Make sure you have set the following environment variables to the external library installation dirs 
* ROOTSYS
* OPENCV_DIR
* BOOST_ROOT
* VLFEAT_INC_DIR, VLFEAT_LIB_DIR: path to VLFEAT headers and libraries 

To build the project:

* Clone this repository into your local $SOURCE_DIR  
  ```git clone https://github.com/simoneriggi/ska-dsh_lmc.git $SOURCE_DIR```
* Enter the source directory $SOURCE_DIR and type:  
  ```make```  
  
