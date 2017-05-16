dnl initial information about the project
AC_INIT([CCJ],[1.0],[jabbari@ualberta.ca])

dnl check if the source foalder is correct
AC_CONFIG_SRCDIR([./CCJ.cpp])

dnl check for C++ preprocessor and compiler and the library compiler
AC_PROG_CXXCPP
AC_PROG_CXX
AC_PROG_CC

dnl special dir for m4
AC_CONFIG_MACRO_DIR([m4])

dnl automake initialisation and check for minimal automake API version 1.9
AM_INIT_AUTOMAKE([1.9])


dnl distribute additional compiler and linker flags
dnl --> set these variables instead of CXXFLAGS or LDFLAGS
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])

## check pkg-config
AC_CHECK_PROG([PKG_CONFIG_AVAILABLE],pkg-config,[true],[false])

dnl get path to simfold
AC_ARG_VAR([SIMFOLD_HOME],
    AC_HELP_STRING(
        [SIMFOLD_HOME],
        [Installation directory of the Simfold library]
      )
)

AS_IF([test "$SIMFOLD_HOME" != "" ],[
    export PKG_CONFIG_PATH="$SIMFOLD_HOME/lib/pkgconfig:$PKG_CONFIG_PATH"
])

dnl configure simfold library via pkg-config
SIMFOLD_OK=false
PKG_CHECK_MODULES([SIMFOLD],[simfold >= 1.0],[SIMFOLD_OK=true])
LIBPATH=${SIMFOLD_LIBS:2:-11}

LIBS="$SIMFOLD_LIBS $LIBS"
AM_CXXFLAGS="$SIMFOLD_CFLAGS $AM_CXXFLAGS"

AS_IF([test "$SIMFOLD_HOME" != "" ],[
    LDFLAGS="-Wl,-rpath,$SIMFOLD_HOME/lib"
])

dnl use the C++ compiler for the following checks
AC_LANG_CPLUSPLUS

dnl TURN ON C++11
AX_CXX_COMPILE_STDCXX_11

dnl Checks for header files.
AC_HEADER_STDC

dnl Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

AS_IF([test "$SIMFOLD_OK" = "true"],[],[
    AC_MSG_ERROR([Cannot find the pkgconfig file of the simfold library.
Try settin the environment variable SIMFOLD_HOME
to the installation path of the simfold library.],
    [-1])
])

AS_IF([test "$SIMFOLD_HOME" = "" ],[
    AC_MSG_WARN([Please set the environment variable SIMFOLD_HOME
to the installation directory of the simfold library.])
])

dnl files to generate via autotools (.am or .in source files)
dnl AC_CONFIG_FILES([src/Makefile])

dnl generate the final Makefile etc.
AC_OUTPUT([Makefile])
