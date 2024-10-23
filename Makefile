# Not a great Makefile.  Specific to my install configuration.  Directories below are easy enough to change.
# The compiler options are taken from the Libmesh examples/eigenproblems Makefile.
LIBMESH_DIR := /scratch/acooksy/libmesh-1.7.1
PETSC_DIR := /scratch/acooksy/petsc-3.13.6
SLEPC_DIR := /scratch/acooksy/slepc-3.13.4
PETSC_ARCH := arch-linux-c-debug
CXX = mpicxx
CFLAGS = -DHAVE_CONFIG_H -std=gnu++17 -DNDEBUG -pthread -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization    -fopenmp 

all: gmat femvib clean

gmat: src/gmat.cpp
	@echo "gmat compilation"
	g++ -Iinclude/gmat/  src/gmat.cpp -o tools/gmat

tune: src/gmat.cpp
	@echo "tune compilation"
	g++ -Iinclude/tune/ src/tune.c -o tools/tune 

femvib: femvib.o
	@echo "FEMvib linking"
	/bin/bash ${LIBMESH_DIR}/libtool  --tag=CXX   --mode=link mpicxx -std=gnu++17 -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization    -fopenmp      -o femvib femvib.o   ${LIBMESH_DIR}/libmesh_opt.la -Wl,-rpath,${SLEPC_DIR}/${PETSC_ARCH}/lib -L${SLEPC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7 -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lslepc -lpetsc -lflapack -lfblas -lpthread -lmpichfort -lmpich -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -ldl ;
# The following attempts to check first if ./deps exists before moving femvib.d there
ifeq (,$(wildcard ./.deps))
	mkdir ./.deps
endif
	mv -f ./femvib.d ./.deps/femvib.Tpo

femvib.o: femvib.cpp
	@echo "FEMvib compilation"
	$(CXX) $(CFLAGS) -I. -I${LIBMESH_DIR}/include. -I/home/acooksy/FEMvib_test/include. -I/scratch/acooksy/eigen-3.4.0  -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I${LIBMESH_DIR}/include -I${LIBMESH_DIR}/contrib/metaphysicl/src/numerics/include -I${LIBMESH_DIR}/contrib/metaphysicl/src/core/include -I${LIBMESH_DIR}/contrib/metaphysicl/src/utilities/include -I${LIBMESH_DIR}/contrib/nanoflann/include -I${LIBMESH_DIR}/contrib/fparser -I${LIBMESH_DIR}/contrib/libHilbert/include -I${LIBMESH_DIR}/contrib/nemesis/v5.22/nemesis -I${LIBMESH_DIR}/contrib/exodusii/v5.22/exodus/cbind/include -I${LIBMESH_DIR}/contrib/netcdf/v4/include -I${LIBMESH_DIR}/contrib/eigen/eigen -I${LIBMESH_DIR}/contrib/gmv -I${LIBMESH_DIR}/contrib/qhull/qhull/src -I${LIBMESH_DIR}/contrib/qhull/qhull/src/libqhullcpp -I${LIBMESH_DIR}/contrib/parmetis/include -I${LIBMESH_DIR}/contrib/metis/include -I${LIBMESH_DIR}/contrib/boost/include -I${LIBMESH_DIR}/contrib/timpi/src/utilities/include -I${LIBMESH_DIR}/contrib/timpi/src/utilities/include -I${LIBMESH_DIR}/contrib/timpi/src/parallel/include -I${LIBMESH_DIR}/contrib/timpi/src/algorithms/include   -MT femvib.o -MD -MP -c -o femvib.o `test -f 'femvib.cpp' || echo './'`femvib.cpp

test:
	if [ -e ./results ] ; then mv ./results ./results_orig ; fi
	@echo "SHO_2D example"
	cp -pr ./examples/results_SHO_2D ./results
	./femvib_run.sh 2
	@echo "SHO 2D eigenvalues"
	head results/eigenvalues.txt
	rm -rf ./results/
	@echo "C6H6O2H2 2D PBC example"
	cp -pr ./examples/results_C6H6O2H2_2D_PBC ./results
	./femvib_run.sh 2
	@echo "C6H6O2H2 2D PBC eigenvalues"
	head results/eigenvalues.txt
	rm -rf ./results/
	@echo "CO2 2D example"
	cp -pr ./examples/results_CO2_2D_str ./results
	./femvib_run.sh 2
	@echo "CO2 2D eigenvalues"
	head results/eigenvalues.txt
	rm -rf ./results/

clean:
	rm -f *.o


