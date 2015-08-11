
###############################################################################
# C++ compiler to use
###############################################################################
#CC = g++
CC = mpic++ -std=c++11 -fpermissive -DUSE_MPI
#CC = c++ -std=c++11
#CC = c++ 
#CC = xlC
#CC = mpCC -DDefineStaticVariables
#CC = KCC
#CC = CC 
#CC = pgCC
##########
# C++ flags
##########
#WFLAGS = -Wall -Werror -ansi -pedantic -v --save-temps 
#OFLAGS = -O3 -funroll-loops 
##########
# Debugging flags to use in compiling
# The -DNDEBUG flag is very important for making the toolset fast.
##########
DEBUGFLAG = -DNDEBUG
#DEBUGFLAG = 
#DEBUGFLAG = -g -DNDEBUG
#DEBUGFLAG = -pg -DNDEBUG

###############################################################################
# FORTRAN compiler to use
###############################################################################
G77 = g77
##########
# FORTRAN flags
##########
G77_SPECIFIC =
G77_FLAGS    =
G77_INCLUDE  = 

##########
# MPI
##########
#MPI_D       = -DUSE_MPI
#MPI_INCLUDE = -I/sw/include
#MPI_LIB     = -L/sw/lib -lmpi -llam

##########
# BLAS
##########
BLAS_D        = -DUSE_BLAS
#BLAS_INCLUDE  = -L/usr/local/i686-Linux/include
#BLAS_LIB      = -L/usr/i686-Linux/lib -llapack -lblas -lg2c
BLAS_LIB      = -llapack -lblas
#BLAS_LIB      = -lcyglapack-0.dll -lcygblas-0.dll

##########
# FFTW
##########
#FFTW_D       = -DUSE_FFTW
#FFTW_INCLUDE = -I/usr/local/i686-Linux/include
#FFTW_LIB     = -L/usr/local/i686-Linux/lib -lfftw3

###############################################################################
# Definitions that rarely need to be changed
###############################################################################

CCFLAGS      = $(WFLAGS) $(DEBUGFLAG) $(OFLAGS) 
INCLUDE      = -I../include \
		$(MPI_D) $(MPI_INCLUDE) \
		$(BLAS_D) $(BLAS_INCLUDE) \
	 	$(FFTW_D) $(FFTW_INCLUDE) 
LIBFLAGS     = $(FFTW_LIB) $(MPI_LIB) $(BLAS_LIB) 
#LIBFLAGS     = -pg

# Default rule, it must come first
all: WangLandau.exe

# General purpose implicit suffix rules
#============================================================
.SUFFIXES:

%.u: %.cpp
	echo $(INCLUDE)
	@echo ------------------------------ Making C++ Dependancy File $@
	$(SHELL) -ec '$(CC) -c $(MAKE_DEPEND_FLAG) $(CC_SPECIFIC) $(CCFLAGS) $(INCLUDE) $(CPPFLAGS) $<; rm -f $*.o'
	-mv -f $*.d $@ 

%.o: %.cpp
	@echo ------------------------------ Making C++ Object File $@
	$(CC) $(CC_SPECIFIC) $(CCFLAGS) $(INCLUDE) -c $< -o $@
	@echo

%.o: %.f
	@echo ------------------------------ Making Fortran Object File $@
	$(G&&) $(G77_SPECIFIC) $(G77_FLAGS) $(G77_INCLUDE) -c $<

%.x: %.o
	@echo ------------------------------ Making X File $@
	$(CC) -o $@ $< $(LIBFLAGS)

%.exe: %.o 
	@echo ------------------------------ Making executable File $@
	$(CC) -o $@ $< $(LIBFLAGS)

hpc:
	mpic++ -DHPC -DUSE_MPI -o Latgas2.exe Latgas2.cpp

# House-cleaning rules
#==================================================
.PHONY : clean cleanObjects
clean:
	-\rm -f *~ *.o *.u *.d

veryclean:
	-\rm -f *~ *.o *.u *.d *.x *.exe

tar:
	cd ../; tar -zcvf caledonia/SRC.tar.gz caledonia/Makefile caledonia/*.cpp caledonia/*.hpp caledonia/fftw3.h
