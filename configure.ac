AC_REVISION([m4_esyscmd_s([git describe --always])])dnl
dnl                                               -*- Autoconf -*-
dnl Process this file with autoreconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([gcrm-io-pnetcdf], [1.0.0], [wkliao@eecs.northwestern.edu])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADER([config.h])
AC_CONFIG_AUX_DIR([./scripts])

dnl Note getting command line should be done before calling AM_INIT_AUTOMAKE
dnl as AM_INIT_AUTOMAKE modifies command line $*
CONFIGURE_ARGS_CLEAN=`echo $* | tr '"' ' '`

dnl parse the version numbers to 3 env variables
VERSION_MAJOR=`echo ${PACKAGE_VERSION} | cut -d. -f1`
VERSION_MINOR=`echo ${PACKAGE_VERSION} | cut -d. -f2`
VERSION_SUB=`echo ${PACKAGE_VERSION} | cut -d. -f3`

dnl Note major, minor, and sub are required
VERSION=$PACKAGE_VERSION

AC_SUBST(VERSION_MAJOR)
AC_SUBST(VERSION_MINOR)
AC_SUBST(VERSION_SUB)
AC_SUBST(VERSION)
AC_SUBST(PACKAGE_VERSION)
AC_SUBST(CONFIGURE_ARGS_CLEAN)

AM_INIT_AUTOMAKE

dnl enable silent rules by default
AM_SILENT_RULES([yes])

dnl configure, Makefile.ins, etc.
AM_MAINTAINER_MODE([enable])

dnl user defined macro for printing messages for debugging
_DEBUG=no
AC_DEFUN([UD_MSG_DEBUG],
         [if test "x${_DEBUG}" = xyes ; then
             AC_MSG_NOTICE(DEBUG: $1)
          fi
         ]
)

AC_ARG_VAR(MPICC, [MPI C compiler, @<:@default: CC@:>@])
ac_user_MPICC=$MPICC
if test "x$MPICC" = x -a "x$CC" != x ; then ac_user_MPICC=$CC ; fi

dnl if MPICC is still empty, search under user's PATH
CANDIDATE_MPICC="${MPICC} mpicc mpicc_r"
CANDIDATE_MPICC+=" mpigcc mpgcc mpigcc_r mpgcc_r"
CANDIDATE_MPICC+=" mpcc_r mpcc mpixlc_r mpixlc"
CANDIDATE_MPICC+=" blrts_xlc mpxlc_r mpxlc"
CANDIDATE_MPICC+=" mpifccpx"
CANDIDATE_MPICC+=" cc"
CANDIDATE_MPICC+=" mpiicc icc"
CANDIDATE_MPICC+=" mpipgcc mppgcc"

dnl check the availability of a list of MPI executables
AC_DEFUN([UD_MPI_PATH_PROGS], [
   ac_mpi_prog_$1=
   AC_PATH_PROGS([ac_mpi_prog_$1], [$2])
   if test "x${ac_mpi_prog_$1}" = x ; then
      AC_PATH_PROGS([ac_mpi_prog_$1], [$2])
   fi
   $1=${ac_mpi_prog_$1}
])

dnl check the availability of one MPI executable in $2
dnl Note $2 can be a compiler name followed by compile options. In this case
dnl we check the first string token, the compiler name.
dnl In addition, $2 can contain the full path of the compiler.
AC_DEFUN([UD_MPI_PATH_PROG], [
   if test "x$2" = x ; then
      AC_MSG_ERROR("2nd argument cannot be NULL")
   else
      AC_MSG_CHECKING($2)
   fi

   dnl 1st token in $2 must be the program name, rests are command-line options
   ac_first_token=`echo $2 | cut -d" " -f1`
   ac_rest_tokens=`echo $2 | cut -d" " -s -f2-`

   ac_mpi_prog_$1=
   AC_PATH_PROG([ac_mpi_prog_$1], [$ac_first_token])

   dnl In case ac_first_token is a full path, the above test may still set
   dnl ac_mpi_prog_$1 to NULL
   if test "x${ac_mpi_prog_$1}" = x ; then
      dnl Note we cannot use AC_CHECK_FILE because it fails for cross compiling
      dnl with error: cannot check for file existence when cross compiling
      ac_mpi_prog_path=`dirname $ac_first_token`
      ac_mpi_prog_name=`basename $ac_first_token`
      if (! test -d "${ac_mpi_prog_path}") ; then
         AC_MSG_ERROR(Directory '${ac_mpi_prog_path}' does not exist)
      fi

      AC_PATH_PROG([ac_mpi_prog_$1], [$ac_mpi_prog_name], [], [$ac_mpi_prog_path])
      if test "x$ac_mpi_prog_$1" = x ; then
         AC_MSG_ERROR($ac_mpi_prog_name cannot be found under $ac_mpi_prog_path)
      fi
   fi
   dnl add back the compile options if there is any
   if test "x$ac_rest_tokens" != x ; then
      ac_mpi_prog_$1="$ac_mpi_prog_$1 $ac_rest_tokens"
   fi

   $1=${ac_mpi_prog_$1}
])

dnl find the full path of MPICC from CANDIDATE_MPICC and MPI_INSTALL
if test "x${ac_user_MPICC}" = x ; then
   dnl if MPICC or CC has not been set by users, then search from
   dnl CANDIDATE_MPICC, and find the full path of MPICC
   UD_MPI_PATH_PROGS([MPICC], [$CANDIDATE_MPICC])
else
   dnl check whether user specified MPICC is valid
   UD_MPI_PATH_PROG([MPICC], [$ac_user_MPICC])
fi

if test "x${MPICC}" = x ; then
   if test "x$ac_user_MPICC" = x ; then
      ERR_MSG="No MPI C compiler can be found"
   else
      ERR_MSG="Specified MPI C compiler \"$ac_user_MPICC\" cannot be found"
   fi
   if test "x$MPI_INSTALL" != x ; then
      ERR_MSG+=" under $MPI_INSTALL"
   fi
   AC_MSG_ERROR([
   -----------------------------------------------------------------------
     $ERR_MSG
     GCRM-IO requires an MPI C compiler. Please specify the location of
     one with the MPICC environment variable.  Abort.
   -----------------------------------------------------------------------])
fi

CC=${MPICC}
AC_PROG_CC
AM_PROG_CC_C_O

dnl enable large file support
AC_SYS_LARGEFILE

dnl check if MPICC works for basic MPI call: MPI_Comm_rank()
AC_CHECK_FUNC([MPI_Comm_rank], [],
   [AC_MSG_ERROR([
    ------------------------------------------------------------
     Invalid MPI compiler specified or detected: "${MPICC}"
     A working MPI compiler is required. Please specify the location
     of one with the MPICC environment variable
    ------------------------------------------------------------])
])

AC_CHECK_FUNC([MPI_File_open], [],
   dnl maybe -lmpio is needed at link stage
   [AC_SEARCH_LIBS([MPI_File_open], [mpio], [],
                   [AC_MSG_ERROR([
   -----------------------------------------------------------------------
     The underneath MPI implementation does not support MPI-IO.
     GCRM-IO requires MPI-IO support to work properly. Abort.
   -----------------------------------------------------------------------])
])])

dnl compute canonical system types
AC_CANONICAL_BUILD
AC_CANONICAL_HOST

PNETCDF_INSTALL=""
AC_ARG_WITH(pnetcdf,
    [AS_HELP_STRING([--with-pnetcdf=/path/to/implementation],
                    [installation prefix for PnetCDF implementation])],
    if test "x${withval}" = xyes; then
        AC_MSG_ERROR(--with-pnetcdf must be given a pathname)
    else
        PNETCDF_INSTALL=${withval}
    fi
)

if test "x${PNETCDF_INSTALL}" != x ; then
    AC_MSG_CHECKING(PnetCDF library)
    CPPFLAGS+=" -I${PNETCDF_INSTALL}/include"
    LDFLAGS+=" -L${PNETCDF_INSTALL}/lib"
fi

have_pnetcdf=no
dnl AC_SEARCH_LIBS([ncmpi_create], [pnetcdf], [have_pnetcdf=yes])
AC_CHECK_LIB([pnetcdf], [ncmpi_create], [have_pnetcdf=yes LIBS+=" -lpnetcdf"])

dnl Checks math library and append -lm instead prepend.
AC_CHECK_LIB([m], [atan2], [LIBS+=" -lm"])

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h sys/time.h unistd.h])

if test "x${have_pnetcdf}" = xyes; then
    AC_CHECK_HEADERS([pnetcdf.h], [], [have_pnetcdf=no])
fi

if test "x${have_pnetcdf}" = xno; then
    AC_MSG_ERROR([
    ------------------------------------------------------------
     The PnetCDF library and header file are required to build
     the GCRM I/O kernel. Use option
         --with-pnetcdf=/path/to/implementation
     to specify the location of PnetCDF build. In addition,
     please make sure the MPI C compiler is compatible with the
     one used to build PnetCDF.    Stopping ...
     Check 'config.log' for more information.
    ------------------------------------------------------------])
fi

dnl AC_SUBST(MPICC)

AC_ARG_ENABLE(debug,
    [AS_HELP_STRING([--enable-debug],
                    [Enable debugging mode. @<:@default: no@:>@])],
    [case "${enableval}" in
             yes) debug=true ;;
             no)  debug=false ;;
             *)   AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;; esac],
    [debug=false]
)
AM_CONDITIONAL(DEBUG, test x"$debug" = x"true")

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_HEADER_TIME

dnl cross compile fails with undefined reference to rpl_realloc and rpl_malloc
dnl AC_FUNC_MALLOC
dnl AC_FUNC_REALLOC

dnl Checks for library functions.
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([gettimeofday memmove memset sqrt strcasecmp strchr strerror strncasecmp])

dnl change default prefix
AC_PREFIX_DEFAULT([${PWD}/run])

AC_PROG_INSTALL
AC_PROG_MAKE_SET

dnl
dnl Check how sed command handling in-place option -i and define SED_I
dnl
AC_DEFUN([UD_PROG_SED_I],
[
   AC_REQUIRE([AC_PROG_SED])
   AC_CACHE_CHECK([for sed handling option -i ], ac_cv_SED_I,[
   cat > conftest.sed_i <<EOF
   test str1
EOF
   ac_cv_err=`$SED -i '' -e 's|str1|str2|g' conftest.sed_i 2>&1`
   if test "x$ac_cv_err" = x ; then
      ac_cv_SED_I="$SED -i ''"
   else
      ac_cv_err=`sed -i'' -e 's|str1|str2|g' conftest.sed_i 2>&1`
      if test "x$ac_cv_err" = x ; then
         ac_cv_SED_I="$SED -i''"
      else
         AC_MSG_ERROR("No proper sed -i option found")
      fi
   fi
   AS_UNSET(ac_cv_err)])
   SED_I="$ac_cv_SED_I"
   AC_SUBST(SED_I)
   rm -f conftest.sed_i
])

dnl AC_PROG_SED and AC_PROG_GREP are only available on autoconf 2.60 and later
AC_PROG_SED
dnl check sed command option -i and set SED_I
UD_PROG_SED_I

dnl tell autoconf the locations of all Makefile templates, .in or .am files
dnl AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 run/Makefile])
AC_OUTPUT
echo \
"--------------------------------------------------------------------
   ${PACKAGE_NAME} Version ${PACKAGE_VERSION}

   Install prefix:  ${prefix}
   MPI C compiler:  ${MPICC}
   CFLAGS:          ${CFLAGS}
   CPPFLAGS:        ${CPPFLAGS}
   LDFLAGS:         ${LDFLAGS}
   LIBS:            ${LIBS}

   Package features:
      GCRM I/O kernel: PnetCDF from '${PNETCDF_INSTALL}'

   Now run 'make' to build the executable gcrm_io in ./run
---------------------------------------------------------------------"
