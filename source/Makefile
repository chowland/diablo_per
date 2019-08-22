# This is the makefile for diablo.
# To compile the code, just type make.  Such an approach makes
# recompilation of the code easy, recompiling only as necessary
# to account for recent changes to the code.
#
# As the user, set the following definitions:

#**********   User Defined Variables Below *********
# Fortran 90 complier to use:
COMPILER = ifort
#COMPILER = gfortran

# Any desired options for the compiler (e.g. -O2, -g, etc.)
# ifort options
USEROPTS = -O2 -mcmodel=large -u
# gfortran options
#USEROPTS = -O2 -mcmodel=large -fimplicit-none

# Location where fftw libraries are installed
#LINKDIR = /nfs/st01/hpc-atmos-jrt51/lib
LINKDIR = /usr/local/lib

# Option to compile with the HDF5 library
HDF5 = TRUE
# **********    END of user definitions ************

# Use the parameters to set flags
COMPOPTS = $(USEROPTS)
#LINKOPTS = -L$(LINKDIR) -lrfftw -lfftw			# to run locally on maths PCs
LINKOPTS = -L$(LINKDIR) -ldrfftw -ldfftw		# to run on Cambridge HPC

MPI = mpi.o mpi_th.o
ifeq ($(HDF5),TRUE)
HDF5_o = write_h5.o
COMPILER = h5pfc
endif

MAIN = diablo.f
HEADER = header

diablo: $(MAIN) diablo_io.o periodic.o \
	write_h5.o fft.o fft_th.o user_rhs.o \
	$(MPI) $(HEADER) grid_def
	$(COMPILER) $(COMPOPTS) $(MAIN) -o diablo \
	diablo_io.o periodic.o \
	write_h5.o fft.o fft_th.o user_rhs.o \
	$(MPI) $(LINKOPTS)

diablo_io.o: diablo_io.f header grid_def
	$(COMPILER) $(COMPOPTS) -c diablo_io.f

periodic.o: periodic.f header grid_def
	$(COMPILER) $(COMPOPTS) -c periodic.f

mpi.o: mpi.f header header_mpi grid_def grid_def_mpi
	$(COMPILER) $(COMPOPTS) -c mpi.f

mpi_th.o: mpi_th.f header header_mpi grid_def grid_def_mpi
	$(COMPILER) $(COMPOPTS) -c mpi_th.f

header : header_mpi grid_def grid_def_mpi

write_h5.o : write_h5.f
	$(COMPILER) $(COMPOPTS) -c write_h5.f

user_rhs.o: user_rhs.f header grid_def grid_def_mpi
	$(COMPILER) $(COMPOPTS) -c user_rhs.f
	
fft.o:  fft.f header grid_def
	$(COMPILER) $(COMPOPTS) -c fft.f

fft_th.o:  fft_th.f header grid_def
	$(COMPILER) $(COMPOPTS) -c fft_th.f

clean:
	rm -f *.o fort.* *~ diablo core

# Compiler specific notes:
#
# Compilation with Absoft Linux Fortran 77 appears to be impossible, as it
# cannot handle the INTEGER*8 option required by FFTW.  If someone finds
# a way around this, please let me know.
#
# Compilation with Absoft Linux Fortran 90 is possible, but the option
# -YEXT_NAMES=LCS must be used as one of the link options so the compiler
# can find the lowercase external library function names.
#
# Compilation with Lahey Fortran 95 (lf95) is possible, but there is an
# underscore incompatability with the FFTW libraries, which are compiled
# with g77.  To get around this, you need to go into fft.f and add
# trailing underscores to the name of every fftw function where they
# appear throughout the code.
