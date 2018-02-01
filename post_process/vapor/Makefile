
# COMPILER = /data/oceanus/ed441/lib/hdf5/1.8.12-ifort/bin/h5fc
COMPILER = h5fc
COMPOPTS = -i-dynamic -fPIC -shared-intel -mcmodel=large
#COMPOPTS = -03 -mcmodel=medium -dhdf5
%COMPILER = h5fc
%COMPOPTS = -i-dynamic -fPIC -shared-intel -mcmodel=medium
# LINKOPTS = -L/usr/local/lib -lrfftw -lfftw -lm
LINKOPTS = -L/mhome/damtp/atmos/jrt51/local/lib -lrfftw -lfftw

h5_to_vapor.o: h5_to_vapor.f fft.o read_h5s.o\
	        header grid_def
	$(COMPILER) $(COMPOPTS) h5_to_vapor.f -o h5_to_vapor \
	fft.o read_h5s.o $(LINKOPTS)

fft.o:  fft.f header grid_def
	$(COMPILER) $(COMPOPTS) -c fft.f

read_h5s.o:	read_h5s.f header grid_def
	$(COMPILER) $(COMPOPTS) -c read_h5s.f

clean:
	rm -f *.o fort.* *~ h5_to_vapor core temp.raw vapor.vdf
	rm -r vapor_data/*
