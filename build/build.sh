#!/bin/bash
export LIBMESH_DIR=/opt/libmesh-1.8.1
export PETSC_DIR=/opt/petsc-3.23.4
export SLEPC_DIR=/opt/slepc-3.23.2
export EIGEN_DIR=/opt/eigen-3.4.0
export PETSC_ARCH=arch-linux-c-debug
export COMPILER=mpicxx
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
	echo ${COMPILER} -std=gnu++17 -DHAVE_CONFIG_H -I. -I${LIBMESH_DIR}/include. -I/home/acooksy/FEMvib/include. -I${EIGEN_DIR} -DNDEBUG  -pthread -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I${LIBMESH_DIR}/include -I${LIBMESH_DIR}/contrib/metaphysicl/src/numerics/include -I${LIBMESH_DIR}/contrib/metaphysicl/src/core/include -I${LIBMESH_DIR}/contrib/metaphysicl/src/utilities/include -I${LIBMESH_DIR}/contrib/nanoflann/include -I${LIBMESH_DIR}/contrib/fparser -I${LIBMESH_DIR}/contrib/libHilbert/include -I${LIBMESH_DIR}/contrib/nemesis/v5.22/nemesis -I${LIBMESH_DIR}/contrib/exodusii/v5.22/exodus/cbind/include -I${LIBMESH_DIR}/contrib/netcdf/v4/include -I${LIBMESH_DIR}/contrib/eigen/eigen -I${LIBMESH_DIR}/contrib/gmv -I${LIBMESH_DIR}/contrib/qhull/qhull/src -I${LIBMESH_DIR}/contrib/qhull/qhull/src/libqhullcpp -I${LIBMESH_DIR}/contrib/parmetis/include -I${LIBMESH_DIR}/contrib/metis/include -I${LIBMESH_DIR}/contrib/boost/include -I${LIBMESH_DIR}/contrib/timpi/src/utilities/include -I${LIBMESH_DIR}/contrib/timpi/src/utilities/include -I${LIBMESH_DIR}/contrib/timpi/src/parallel/include -I${LIBMESH_DIR}/contrib/timpi/src/algorithms/include   -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization    -fopenmp    -MT femvib.o -MD -MP -c -o femvib.o 
	# -fpermissive added against compiler advice -- there is deprecated syntax in use of MAX in src/include/krig.h line 71
	${COMPILER} -std=gnu++17 -fpermissive -DHAVE_CONFIG_H -I. -I${LIBMESH_DIR}/include. -I/home/acooksy/FEMvib/include. -I${EIGEN_DIR} -DNDEBUG  -pthread -I${SLEPC_DIR}/include -I${SLEPC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I${LIBMESH_DIR}/include -I${LIBMESH_DIR}/contrib/metaphysicl/src/numerics/include -I${LIBMESH_DIR}/contrib/metaphysicl/src/core/include -I${LIBMESH_DIR}/contrib/metaphysicl/src/utilities/include -I${LIBMESH_DIR}/contrib/nanoflann/include -I${LIBMESH_DIR}/contrib/fparser -I${LIBMESH_DIR}/contrib/libHilbert/include -I${LIBMESH_DIR}/contrib/nemesis/v5.22/nemesis -I${LIBMESH_DIR}/contrib/exodusii/v5.22/exodus/cbind/include -I${LIBMESH_DIR}/contrib/netcdf/v4/include -I${LIBMESH_DIR}/contrib/eigen/eigen -I${LIBMESH_DIR}/contrib/gmv -I${LIBMESH_DIR}/contrib/qhull/qhull/src -I${LIBMESH_DIR}/contrib/qhull/qhull/src/libqhullcpp -I${LIBMESH_DIR}/contrib/parmetis/include -I${LIBMESH_DIR}/contrib/metis/include -I${LIBMESH_DIR}/contrib/boost/include -I${LIBMESH_DIR}/contrib/timpi/src/utilities/include -I${LIBMESH_DIR}/contrib/timpi/src/utilities/include -I${LIBMESH_DIR}/contrib/timpi/src/parallel/include -I${LIBMESH_DIR}/contrib/timpi/src/algorithms/include   -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization    -fopenmp    -MT femvib.o -MD -MP -c -o femvib.o `test -f 'src/femvib.cpp' || echo './'`src/femvib.cpp
	/bin/bash ${LIBMESH_DIR}/libtool  --tag=CXX   --mode=link mpicxx -std=gnu++17 -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization    -fopenmp      -o femvib femvib.o   ${LIBMESH_DIR}/libmesh_opt.la -Wl,-rpath,${SLEPC_DIR}/${PETSC_ARCH}/lib -L${SLEPC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib -L${PETSC_DIR}/${PETSC_ARCH}/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7 -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lslepc -lpetsc -lflapack -lfblas -lpthread -lmpichfort -lmpich -lgfortran -lm -lgcc_s -lquadmath -lstdc++ -ldl ;
	rm -f femvib.o
	if [ ! -e ./.deps ] ; then
		mkdir ./.deps
	fi
	mv ./femvib.d .deps/femvib_ep-ex2.Tpo
fi

# OPTION DEFINITIONS
# -MP, -MT, -MD. -MF are all preprocessor options to gcc;  see https://gcc.gnu.org/onlinedocs/gcc/Preprocessor-Options.html
# -MP  This option instructs CPP to add a phony target for each dependency other than the main file, causing each to depend on nothing. These dummy rules work around errors make gives if you remove header files without updating the Makefile to match.
# -MT [target]  Change the target of the rule emitted by dependency generation. By default CPP takes the name of the main input file, deletes any directory components and any file suffix such as ‘.c’, and appends the platform’s usual object suffix. The result is the target.  An -MT option sets the target to be exactly the string you specify. If you want multiple targets, you can specify them as a single argument to -MT, or use multiple -MT options. 
# -MF [filename]   When used with `-M' or `-MM', specifies a file to write the dependencies to.  If no `-MF' switch is given the preprocessor sends the rules to the same place it would have sent preprocessed output.  When used with the driver options `-MD' or `-MMD', `-MF' overrides the default dependency output file.
# -MD  equivalent to `-M -MF FILE', except that `-E' is not implied.  The driver determines FILE based on whether an `-o' option is given.  If it is, the driver uses its argument but with a suffix of `.d', otherwise it take the basename of the input file and applies a `.d' suffix.  If `-MD' is used in conjunction with `-E', any `-o' switch is understood to specify the dependency output file (but *note -MF: dashMF.), but if used without `-E', each `-o' is understood to specify a target object file.  Since `-E' is not implied, `-MD' can be used to generate a dependency output file as a side-effect of the compilation process.
#
# -DNDEBUG   -D is the GCC option to define a macro. NDEBUG is the macro to be defined to turn off asserts as mandated by the C standard.
#
# -pthread   This causes files to be compiled with -D_REENTRANT, and linked with -lpthread.  Using _REENTRANT, on GNU libc, changes the way some libc headers work. As a specific example, it makes errno call a function returning a thread-local location.
