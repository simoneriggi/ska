#
# Makefile for CAESAR tool
#
#
###  BASE DIRS
TOPDIR= $(shell pwd)
BINDIR  = $(TOPDIR)/bin
LIBDIR  = $(TOPDIR)/lib
SRCDIR  = $(TOPDIR)/src
INCDIR  = $(TOPDIR)/include
OBJDIR  = $(TOPDIR)/obj

### SOURCE FILES
USER_SRCS = $(wildcard $(SRCDIR)/*.cc)
files := $(foreach USER_SRCS,$(SRCDIR),$(USER_SRCS))

### HEADERS USED IN DICTIONARY GENERATION
HEADERS_DICT = $(INCDIR)/ImgFITSReader.h $(INCDIR)/FITSReader.h $(INCDIR)/Img.h $(INCDIR)/Source.h $(INCDIR)/Contour.h $(INCDIR)/Region.h $(INCDIR)/BkgFinder.h $(INCDIR)/ChanVeseSegmentation.h $(INCDIR)/EllipseUtils.h $(INCDIR)/SourceFitter.h $(INCDIR)/SLICSegmenter.h $(INCDIR)/SLICUtils.h $(INCDIR)/HClust.h $(INCDIR)/OutlierDetector.h $(INCDIR)/Interpolator.h $(INCDIR)/BkgFinder.h $(INCDIR)/ZernikeMoments.h

### OBJECT GENERATION
OBJS = $(USER_SRCS:.cc=.o)
MAIN_FINDSOURCE = FindSource.cc
MAIN_FINDSOURCE_OBJ = $(MAIN_FINDSOURCE:.cc=.o)
MAIN_FINDBKG = FindBackground.cc
MAIN_FINDBKG_OBJ = $(MAIN_FINDBKG:.cc=.o)

## Get platform 32/64 bit
LBITS   = $(shell getconf LONG_BIT)

## Get SVN revision number
SVN_REVISION_VERSION = -D'SVN_REVISION_VERSION="$(shell svnversion -n .)"'

# Set executable a name
EXE = libImg.so
FINDSOURCE_EXE = FindSource
FINDBKG_EXE = FindBackground
#
#############################################################

## You should not need to change anything below this line ###

.PHONY: all depend clean


######################################
###  CPPFLAGS & CXXFLAGS  ############
######################################

CPPFLAGS = -I$(INCDIR)

ifeq ($(LBITS),64)
  # do 64 bit stuff here
	CPPFLAGS += -I/usr/include -pthread -m64
	CXXFLAGS = -std=c++11 -O2 -Wall -Wunused -Wuninitialized -fPIC -pthread -m64 
	SYSLIBDIR = /usr/lib/x86_64-linux-gnu
else
  # do 32 bit stuff here
	CPPFLAGS += -I/usr/include -pthread -m32
  CXXFLAGS = -std=c++11 -O2 -Wall -Wunused -Wuninitialized -fPIC -pthread -m32 
	SYSLIBDIR = /usr/lib
endif

SOFLAGS = -fPIC -ggdb3 -Wall -shared

ROOT_CPPFLAGS = -I$(ROOTSYS)/include
BOOST_CPPFLAGS = -I$(BOOSTROOT)/include
OPENCV_CPPFLAGS = $(shell pkg-config $(OPENCV_DIR)/lib/pkgconfig/opencv.pc --cflags)

#VLFEAT_DIR= /home/riggi/Software/VLFeat/lib
#VLFEAT_CPPFLAGS = -I$(VLFEAT_DIR)
VLFEAT_CPPFLAGS = -I$(VLFEAT_SRC_DIR)

R_HOME = $(shell R RHOME)
R_CPPFLAGS = $(shell $(R_HOME)/bin/R CMD config --cppflags)
RCPP_CPPFLAGS = $(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RINSIDE_CPPFLAGS = $(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)

CPPFLAGS += $(ROOT_CPPFLAGS) $(WCSTOOLS_CPPFLAGS) $(OPENCV_CPPFLAGS) $(VLFEAT_CPPFLAGS) $(BOOST_CPPFLAGS)
CPPFLAGS += $(R_CPPFLAGS) $(RCPP_CPPFLAGS) $(RINSIDE_CPPFLAGS) # Add R flags

###########################
###  LDFLAGS   ############
###########################
ROOT_LDFLAGS = $(shell root-config --libs) -L$(ROOTSYS)/lib -lSpectrum
ROOT_LDFLAGS += -lFITSIO
BOOST_LDFLAGS = -L$(BOOSTROOT)/lib -lboost_filesystem -lboost_system
SYSTEM_LDFLAGS = -L$(SYSLIBDIR) -lrt 
OPENCV_LDFLAGS = $(shell pkg-config $(OPENCV_DIR)/lib/pkgconfig/opencv.pc --libs)
##VLFEAT_LDFLAGS = $(VLFEAT_DIR)/libvtfeat.a
VLFEAT_LDFLAGS = -L$(VLFEAT_LIB_DIR) -lvl


R_LDFLAGS = $(shell $(R_HOME)/bin/R CMD config --ldflags)
RBLAS_LDFLAGS = $(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
RLAPACK_LDFLAGS = $(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)
RCPP_LDFLAGS = $(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)
RINSIDE_LDFLAGS = $(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)
GSL_LDFLAGS = -L/usr/lib -lgsl -lblas

LDFLAGS = $(ROOT_LDFLAGS) $(BOOST_LDFLAGS) $(SYSTEM_LDFLAGS) $(WCSTOOLS_LDFLAGS) $(OPENCV_LDFLAGS) $(VLFEAT_LDFLAGS)
LDFLAGS += $(R_LDFLAGS) $(RBLAS_LDFLAGS) $(RLAPACK_LDFLAGS) $(RCPP_LDFLAGS) $(RINSIDE_LDFLAGS) $(GSL_LDFLAGS)

################################################################


all: GETOBJS PRINTINFO $(EXE) PUTOBJS $(FINDSOURCE_EXE) $(FINDBKG_EXE)


PRINTINFO: 
	@echo 'Compiling $(SVN_REVISION_VERSION) trunk version on a $(LBITS) bit machine' \

REMOVEOBJS:
	- rm -f *.o $(SRCDIR)/*.o 

GETOBJS:
	@echo "Put object and lib files again to $(SRCDIR) dir"
	- mv -f $(OBJDIR)/*.o $(SRCDIR) 2>/dev/null
	- mv -f $(LIBDIR)/*.so $(SRCDIR) 2>/dev/null


PUTOBJS:
	@echo "Moving object files to $(OBJDIR) dir and libraries to $(LIB) dir"
	- mv -f $(SRCDIR)/*.o $(OBJDIR) 2>/dev/null
	- mv -f *.o $(OBJDIR) 2>/dev/null
	- mv -f $(EXE) $(LIBDIR) 2>/dev/null

ClassDictionary.cc: $(HEADERS_DICT) LinkDef.h
	rootcint -f $@ -c -p $(CPPFLAGS) $^
#rootcint $(CXXFLAGS) -f $@ -c $^
	
$(EXE): $(OBJS) ClassDictionary.o
	@echo "Compiler is $(CXX) or $(CC), options are $(CXXFLAGS), generating $@ shared library..."
	@$(CXX) $(CXXFLAGS) $(SOFLAGS) $(CPPFLAGS) $^ -o $@ $(LDFLAGS) 


$(FINDSOURCE_EXE): $(MAIN_FINDSOURCE_OBJ)
	@echo "Building find source pipeline ..."
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(MAIN) $(LDFLAGS) -L$(LIBDIR) -lImg

$(FINDBKG_EXE): $(MAIN_FINDBKG_OBJ)
	@echo "Building find background pipeline ..."
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(MAIN) $(LDFLAGS) -L$(LIBDIR) -lImg

#############################################################
# gcc can generate the dependency list

depend: Make-depend

Make-depend: $(USER_SRCS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $^ > $@

clean:
	- rm -f $(EXE) $(FINDSOURCE_EXE) $(FINDBKG_EXE) *.o $(OBJDIR)/*.o $(SRCDIR)/*.o ClassDictionary.cc ClassDictionary.h ClassDictionary_rdict.pcm $(BINDIR)/*.o *.so $(LIBDIR)/*.so *.ps core Make-depend


#############################################################
# 'make run' will run the thing



#############################################################
# the lines below are for running with debugger 'make run_gdb'

.INTERMEDIATE: gdb.cmdl

# batch mode gdb needs a file with commands
gdb.cmdl:
	echo "r -b bootstrap.xml" > $@

run_gdb: gdb.cmdl $(EXE)
	gdb -batch -x $< ./$(EXE) && touch $@

#include Make-depend
