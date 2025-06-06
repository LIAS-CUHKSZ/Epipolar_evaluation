VERSION     =  7.3.18
SDPA_DIR    =  /usr/local
SDPA_LIB    =  /usr/local/lib/libsdpa.a
MUMPS_DIR   =  /home/neo/sdpa-7.3.18/mumps/build
CC          =  gcc
CFLAGS      =  -g -O2 -Wall -fPIC -funroll-all-loops
CXX         =  g++
CXXFLAGS    =  -g -O2 -D_REENTRANT -Wall -fPIC -funroll-all-loops
CPPFLAGS    = 
FC          =  gfortran
FCFLAGS     =  -g -O2 -Wall -fPIC -funroll-all-loops
LAPACK_LIBS =  -llapack
BLAS_LIBS   =  -lblas
PTHREAD_INCLUDE=  -I/usr/include
PTHREAD_LIBS=  -lpthread
FCLIBS      =  -L/usr/lib/gcc/x86_64-linux-gnu/11 -L/usr/lib/gcc/x86_64-linux-gnu/11/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/11/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/11/../../.. -lgfortran -lm -lquadmath
MUMPS_INCLUDE =  -I/home/neo/sdpa-7.3.18/mumps/build/include
MUMPS_LIBS =  -L/home/neo/sdpa-7.3.18/mumps/build/lib -ldmumps -lmumps_common -lpord -L/home/neo/sdpa-7.3.18/mumps/build/libseq -lmpiseq
SDPA_LIBS =  /usr/local/lib/libsdpa.a -L/home/neo/sdpa-7.3.18/mumps/build/lib -ldmumps -lmumps_common -lpord -L/home/neo/sdpa-7.3.18/mumps/build/libseq -lmpiseq -llapack -lblas -L/usr/lib/gcc/x86_64-linux-gnu/11 -L/usr/lib/gcc/x86_64-linux-gnu/11/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/11/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/11/../../.. -lgfortran -lm -lquadmath -lpthread
