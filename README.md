# diablo_per
Triply periodic version of Diablo DNS program (Fortran)

### Prerequisites
1. Fortran compiler (gfortran & ifort tested)
2. MPI implementation (MPICH 3 tested)
3. HDF5 parallel build (1.8.xx tested)
4. FFTW 2 (subroutines need updating before FFTW 3 can be used)

#### Quick and easy setup guide from scratch on Ubuntu
- Install gfortran and mpich using apt-get
```
sudo apt-get install gfortran
sudo apt-get install mpich
```
- Download HDF5 1.8.xx source code & build using mpich wrapper for gfortran
```
cd hdf5-1.8.xx
CC=/usr/bin/mpicc
FC=/usr/bin/mpif90
./configure --enable-parallel --prefix=/usr/local/hdf5 --enable-fortran
make
sudo make install
```
- Add the following line to your .bashrc file in the home directory
```
PATH=/usr/local/hdf5/bin:$PATH ; export PATH
```
- Download FFTW 2 source code & build libraries using hdf5 wrapper for gfortran
```
cd fftw-2.1.5
./configure F77=h5pfc
make
sudo make install
```
