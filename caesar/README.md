# CAESAR

Compact And Extended Source Automated Recognition

##**About**  
CAESAR is a C++ software tool for automated source finding in astronomical maps. 
It is distributed under the GNU General Public License v3.0. 
If you use CAESAR for your research, it would be nice to acknowledge it in your papers by citing the following paper:

* S. Riggi et al., "Automated detection of extended sources in radio maps:
progress from the SCORPIO survey", MNRAS (2016) doi: 10.1093/mnras/stw982, arXiv:1605.01852

##**Status**
At present CAESAR is still being developed with major functionalities to be integrated in this repository. The API could be therefore subjected to significant changes. 

##**Installation**  
Install the project dependencies:  
* ROOT [https://root.cern.ch/], to be built with FITSIO, PyROOT, RInterface options enabled
* OpenCV [http://opencv.org/]
* R [https://www.r-project.org/], install also these additional packages: RInside, Rcpp, rrcovHD, truncnorm, FNN, akima
* log4cxx [https://logging.apache.org/log4cxx/]
* boost [http://www.boost.org/]
* python (>=2.7) [https://www.python.org/], install also these additional modules: pyfits
* cmake (>=2.8) [https://cmake.org]  
  
Make sure you have set the following environment variables to the external library installation dirs 
* ROOTSYS
* OPENCV_DIR
* BOOST_ROOT
* LOG4CXX_ROOT

Add also the following paths to the PKG_CONFIG_PATH environment var: 
* $LOG4CXX_ROOT/lib/pkgconfig  

CAESAR depends also on the wcstools and linterp libraries which are already provided in the external/ directory. Note that the provided wcslib was slightly modified with respect to the original release to avoid naming conflicts with the R package.

cmake should find all needed include dirs and libraries used to build the project.

To build and install the project:

* Clone this repository into your local $SOURCE_DIR  
  ```git clone https://github.com/simoneriggi/ska-dsh_lmc.git $SOURCE_DIR```
* Create the build and install directories: $BUILD_DIR, $INSTALL_DIR  
* In the build directory:  
  ```cmake -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR $SOURCE_DIR```  
  ```make```  
  ```make install```  
  
