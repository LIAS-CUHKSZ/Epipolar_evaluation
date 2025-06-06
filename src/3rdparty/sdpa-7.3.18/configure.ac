dnl
dnl Configure for SDPA 
dnl Copyright (C) 2004-2024 SDPA Project
dnl

AC_INIT([sdpa],[7.3.18])
AM_INIT_AUTOMAKE
AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_CXX
AC_PROG_FC
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS
AC_CANONICAL_HOST

AC_ARG_WITH(pthread-include,
[  --with-pthread-include  Gives pthread include you use(e.g., "-I/usr/include")],
PTHREAD_INCLUDE=$withval
)

AC_ARG_WITH(pthread-libs,
[  --with-pthread-libs     Gives pthread name you use(e.g., "-lpthread")],
PTHREAD_LIBS=$withval
)

AC_ARG_WITH(blas,
[  --with-blas             Gives BLAS full name you use(e.g., "/usr/lib/liblblas.a")],
BLAS_LIBS=$withval
)

AC_ARG_WITH(lapack,
[  --with-lapack           Gives LAPACK full name you use(e.g., "/usr/lib/liblapack.a")],
LAPACK_LIBS=$withval
)

AC_ARG_WITH(mumps-include,
[  --with-mumps-include    Gives MUMPS include you use(e.g., "-I/usr/include    or    -I/usr/lib/mumps/include")],
MUMPS_INCLUDE=$withval
)

AC_ARG_WITH(mumps-libs,
[  --with-mumps-libs       Gives MUMPS library you use(e.g., "/usr/lib/libdmumps_seq.a /usr/lib/libmumps_common_seq.a /usr/lib/libmpiseq_seq.a /usr/lib/libpord_seq.a    or   /usr/lib/mumps/lib/libdmumps.a /usr/lib/mumps/lib/libmumps_common.a /usr/lib/mumps/lib/libpord.a /usr/lib/mumps/libseq/libmpiseq.a")],
MUMPS_LIBS=$withval
)

dnl spooles requires -lm when linking.
LIB="-lm"

dnl How thread we link thread library esp. needed for spooles.
case "$host" in
        *-*-linux*|*-*-uclinux*)
            pthread_cflags="-D_REENTRANT"
            ;;
esac
CXXFLAGS="$CXXFLAGS $pthread_cflags"

dnl -fPIC is necessary to link with mex files on 64bit linux
CFLAGS="$CFLAGS -Wall -fPIC -funroll-all-loops"
CXXFLAGS="$CXXFLAGS -Wall -fPIC -funroll-all-loops"
FCFLAGS="$FCFLAGS -Wall -fPIC -funroll-all-loops"

AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(CPPFLAGS)
AC_SUBST(FCFLAGS)
AC_SUBST(FCLIBS)


dnl check pthread rountines
HAVE_PTHREAD=""
dnl default blas name is libpthread, thus we pass -lpthread.
if test x"$PTHREAD_INCLUDE" = x; then
  PTHREAD_INCLUDE="-I/usr/include"
fi
if test x"$PTHREAD_LIBS" = x; then
  PTHREAD_LIBS="-lpthread"
fi
save_LIBS="$LIBS"
LIBS="$LIBS $PTHREAD_LIBS"
checkLIBS=$LIBS
AC_CHECK_FUNC(pthread_create, [HAVE_PTHREAD=yes])
LIBS="$save_LIBS"
AC_SUBST(PTHREAD_INCLUDE)
AC_SUBST(PTHREAD_LIBS)

if test x"$HAVE_PTHREAD" != xyes; then
  echo "Pthread function (pthread_create) was compiled by"
  echo "( $CC $CXXFLAGS $PTHREAD_INCLUDE conftest.c $checkLIBS )"
  echo "but failed."
  echo "(1) Please verify config.log to check whether libraries are linked correctrly"
  AC_MSG_ERROR([Pthread link failed])
fi


dnl check BLAS/LAPACK rountines
HAVE_BLAS=""

dnl default blas name is libblas, thus we pass -lblas.
if test x"$BLAS_LIBS" = x; then
  BLAS_LIBS="-lblas"
fi
save_LIBS="$LIBS"
AC_FC_FUNC(dgemm)
LIBS="$LIBS $BLAS_LIBS $FCLIBS $PTHREAD_LIBS"
checkLIBS=$LIBS
AC_CHECK_FUNC($dgemm, [HAVE_BLAS=yes])
LIBS="$save_LIBS"
AC_SUBST(BLAS_LIBS)

if test x"$HAVE_BLAS" != xyes; then
  echo "BLAS function (dgemm) was compiled by"
  echo "( $CC $CFLAGS conftest.c $checkLIBS )"
  echo "but failed."
  echo "(1) Please verify config.log to check whether libraries are linked correctrly"
  echo "(2) Please install one of your favorite BLAS"
  echo "OpenBLAS:  http://xianyi.github.com/OpenBLAS/"
  echo "BLAS:      http://www.netlib.org/blas/"
  echo "ATLAS:     http://math-atlas.sourceforge.net/"
  AC_MSG_ERROR([BLAS link failed])
fi
          
dnl Check whether we can call LAPACK subroutines
HAVE_LAPACK=""

dnl default lapack name is liblapack, thus we pass -llapack.
if test x"$LAPACK_LIBS" = x; then
  LAPACK_LIBS="-llapack"
fi
save_LIBS="$LIBS"
AC_FC_FUNC(dsyev)
LIBS="$LIBS $LAPACK_LIBS $BLAS_LIBS $FCLIBS $PTHREAD_LIBS"
checkLIBS=$LIBS
AC_CHECK_FUNC($dsyev, [HAVE_LAPACK=yes])
LIBS="$save_LIBS"
AC_SUBST(LAPACK_LIBS)

if test x"$HAVE_LAPACK" != xyes; then
  echo "LAPACK function (dsyev) was compiled by"
  echo "( $CC $CFLAGS conftest.c $checkLIBS )"
  echo "but failed."
  echo "Please verify config.log to check whether libraries are linked correctrly"
  echo "Please install LAPACK"
  echo "LAPACK:   http://www.netlib.org/lapack/"
  AC_MSG_ERROR([LAPACK link failed])
fi

dnl Check whether we can call MUMPS subroutines
HAVE_MUMPS=""
dnl default lapack name is liblapack, thus we pass -llapack.
if test x"$MUMPS_INCLUDE" = x; then
  MUMPS_INCLUDE="-I/usr/include"
fi
if test x"$MUMPS_LIBS" = x; then
   MUMPS_LIBS="/usr/lib/libdmumps_seq.a /usr/lib/libmumps_common_seq.a /usr/lib/libmpiseq_seq.a /usr/lib/libpord_seq.a"
fi
save_LIBS="$LIBS"
AC_FC_FUNC(dmumps)
LIBS="$LIBS $MUMPS_LIBS $LAPACK_LIBS $BLAS_LIBS $FCLIBS $PTHREAD_LIBS"
checkLIBS=$LIBS
AC_CHECK_FUNC($dmumps, [HAVE_MUMPS=yes])
LIBS="$save_LIBS"

if test x"$HAVE_MUMPS" != xyes; then
  echo "MUMPS link failed"
  echo "MUMPS function (dmumps) was compiled by"
  echo "( $CC $CFLAGS $MUMPS_INCLUDE conftest.c $checkLIBS )"
  echo "but failed."
  echo "MUMPS packages will be downloaded by wget command."
  MUMPS_DIR=`pwd`/mumps/build
  MUMPS_INCLUDE="-I${MUMPS_DIR}/include"
  MUMPS_LIBS="-L${MUMPS_DIR}/lib -ldmumps -lmumps_common -lpord -L${MUMPS_DIR}/libseq -lmpiseq"
fi
AC_SUBST(MUMPS_INCLUDE)
AC_SUBST(MUMPS_LIBS)
AC_SUBST(MUMPS_DIR)
AM_CONDITIONAL(HAVE_PACKAGE_MUMPS, [test x"$HAVE_MUMPS" = xyes])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
