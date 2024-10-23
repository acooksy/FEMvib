#!/bin/bash
export LIBMESH_DIR=/scratch/acooksy/libmesh-1.7.1
export PETSC_DIR=/scratch/acooksy/petsc-3.13.6
export SLEPC_DIR=/scratch/acooksy/slepc-3.13.4
export PETSC_ARCH=arch-linux-c-debug
rootdir=`pwd`
cd $rootdir
if [ ! -e ./tools/gmat ] ; then
	echo "compiling gmat" 
	g++ -Iinclude/gmat/ -o tools/gmat src/gmat.cpp ;
fi
if [ ! -e ./tools/tune ] ; then
	echo "compiling tune" 
	g++ -I./include/tune/ -o tools/tune src/tune.c
fi
if [ ! -e ./femvib ] ; then
	echo "compiling femvib" 
#	mpicxx -std=gnu++17 -DHAVE_CONFIG_H -I. -I/scratch/acooksy/libmesh-1.7.1/include. -I/home/acooksy/FEMvib_test/include. -I/scratch/acooksy/eigen-3.4.0  -DNDEBUG  -pthread -I/scratch/acooksy/slepc-3.13.4/include -I/scratch/acooksy/slepc-3.13.4/arch-linux-c-debug/include -I/scratch/acooksy/petsc-3.13.6/include -I/scratch/acooksy/petsc-3.13.6/arch-linux-c-debug/include -I/scratch/acooksy/libmesh-1.7.1/include -I/scratch/acooksy/libmesh-1.7.1/contrib/metaphysicl/src/numerics/include -I/scratch/acooksy/libmesh-1.7.1/contrib/metaphysicl/src/core/include -I/scratch/acooksy/libmesh-1.7.1/contrib/metaphysicl/src/utilities/include -I/scratch/acooksy/libmesh-1.7.1/contrib/nanoflann/include -I/scratch/acooksy/libmesh-1.7.1/contrib/fparser -I/scratch/acooksy/libmesh-1.7.1/contrib/libHilbert/include -I/scratch/acooksy/libmesh-1.7.1/contrib/nemesis/v5.22/nemesis -I/scratch/acooksy/libmesh-1.7.1/contrib/exodusii/v5.22/exodus/cbind/include -I/scratch/acooksy/libmesh-1.7.1/contrib/netcdf/v4/include -I/scratch/acooksy/libmesh-1.7.1/contrib/eigen/eigen -I/scratch/acooksy/libmesh-1.7.1/contrib/gmv -I/scratch/acooksy/libmesh-1.7.1/contrib/qhull/qhull/src -I/scratch/acooksy/libmesh-1.7.1/contrib/qhull/qhull/src/libqhullcpp -I/scratch/acooksy/libmesh-1.7.1/contrib/parmetis/include -I/scratch/acooksy/libmesh-1.7.1/contrib/metis/include -I/scratch/acooksy/libmesh-1.7.1/contrib/boost/include -I/scratch/acooksy/libmesh-1.7.1/contrib/timpi/src/utilities/include -I/scratch/acooksy/libmesh-1.7.1/contrib/timpi/src/utilities/include -I/scratch/acooksy/libmesh-1.7.1/contrib/timpi/src/parallel/include -I/scratch/acooksy/libmesh-1.7.1/contrib/timpi/src/algorithms/include   -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization    -fopenmp    -MT femvib.o -MD -MP -c -o femvib.o `test -f 'femvib.C' || echo './'`femvib.C
	mpicxx -std=gnu++17 -DHAVE_CONFIG_H -I. -I${LIBMESH_DIR}/include. -I/home/acooksy/FEMvib_test/include. -I/scratch/acooksy/eigen-3.4.0  -DNDEBUG  -pthread -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I${LIBMESH_DIR}/include -I${LIBMESH_DIR}/contrib/metaphysicl/src/numerics/include -I${LIBMESH_DIR}/contrib/metaphysicl/src/core/include -I${LIBMESH_DIR}/contrib/metaphysicl/src/utilities/include -I${LIBMESH_DIR}/contrib/nanoflann/include -I${LIBMESH_DIR}/contrib/fparser -I${LIBMESH_DIR}/contrib/libHilbert/include -I${LIBMESH_DIR}/contrib/nemesis/v5.22/nemesis -I${LIBMESH_DIR}/contrib/exodusii/v5.22/exodus/cbind/include -I${LIBMESH_DIR}/contrib/netcdf/v4/include -I${LIBMESH_DIR}/contrib/eigen/eigen -I${LIBMESH_DIR}/contrib/gmv -I${LIBMESH_DIR}/contrib/qhull/qhull/src -I${LIBMESH_DIR}/contrib/qhull/qhull/src/libqhullcpp -I${LIBMESH_DIR}/contrib/parmetis/include -I${LIBMESH_DIR}/contrib/metis/include -I${LIBMESH_DIR}/contrib/boost/include -I${LIBMESH_DIR}/contrib/timpi/src/utilities/include -I${LIBMESH_DIR}/contrib/timpi/src/utilities/include -I${LIBMESH_DIR}/contrib/timpi/src/parallel/include -I${LIBMESH_DIR}/contrib/timpi/src/algorithms/include   -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization    -fopenmp    -MT femvib.o -MD -MP -c -o femvib.o `test -f 'femvib.C' || echo './'`femvib.C
#	/bin/bash /scratch/acooksy/libmesh-1.7.1/libtool  --tag=CXX   --mode=link mpicxx -std=gnu++17 -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization    -fopenmp      -o femvib femvib.o   /scratch/acooksy/libmesh-1.7.1/libmesh_opt.la -Wl,-rpath,/scratch/acooksy/slepc-3.13.4/arch-linux-c-debug/lib -L/scratch/acooksy/slepc-3.13.4/arch-linux-c-debug/lib -Wl,-rpath,/scratch/acooksy/petsc-3.13.6/arch-linux-c-debug/lib -L/scratch/acooksy/petsc-3.13.6/arch-linux-c-debug/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7 -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lslepc -lpetsc -lflapack -lfblas -lpthread -lmpichfort -lmpich -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -ldl ;
	/bin/bash ${LIBMESH_DIR}/libtool  --tag=CXX   --mode=link mpicxx -std=gnu++17 -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization    -fopenmp      -o femvib femvib.o   ${LIBMESH_DIR}/libmesh_opt.la -Wl,-rpath,${SLEPC_DIR}/${PETSC_ARCH}/lib -L${SLEPC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7 -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lslepc -lpetsc -lflapack -lfblas -lpthread -lmpichfort -lmpich -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -ldl ;
	rm -f femvib.o
	if [ ! -e ./.deps ] ; then
		mkdir ./.deps
	fi
	mv ./femvib.d .deps/femvib_ep-ex2.Tpo
fi
