# This is an example makefile for the user defined refinement in HMR
# Matt, please convert this into a Makefile.in and figure out how
# to provide it to the user
TPLS=/home/maute/tpls
MORISROOT=/home/messe/codes/moris
CXX=$(TPLS)/openmpi/gcc/bin/mpic++

BOOST   = -I$(TPLS)/boost/include
LAPACK  = -DF77ADD_ -DMORIS_HAVE_MKL -I$(TPLS)/mkl/include
LINALG  = -DMORIS_USE_ARMA -I$(TPLS)/armadillo/include
#LINALG  = -DMORIS_USE_EIGEN  -I$(TPLS)/eigen/include

MORIS = -I$(MORISROOT)/projects \
	-I$(MORISROOT)/snippets \
	-I$(MORISROOT)/include \
	-I$(MORISROOT)/projects/MRS/ASR/src \
	-I$(MORISROOT)/projects/MRS/CHR/src \
	-I$(MORISROOT)/projects/MRS/CNT/src \
	-I$(MORISROOT)/projects/MRS/COR/src \
	-I$(MORISROOT)/projects/MRS/EXC/src \
	-I$(MORISROOT)/projects/MRS/IOS/src \
	-I$(MORISROOT)/projects/LINALG/src \
	-I$(MORISROOT)/projects/COM/src \
	-I$(MORISROOT)/projects/MTK/src \
	-I$(MORISROOT)/projects/HMR/src
 
INCLUDE = $(BOOST) $(LAPACK) $(LINALG) $(MORIS)         

DFLAGS   = -DMORIS_USE_32BIT
CXXFLAGS = -O2 -fPIC -m64 -std=c++11 -g

hmr_userdef.so: /home/messe/codes/moris/projects/HMR/tutorials/HMR_User_Defined_Refinement.cpp
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(INCLUDE) -shared $< -o hmr_userdef.so

clean:
	rm -v hmr_userdef.so
