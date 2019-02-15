cmake_minimum_required(VERSION 3.0) # setting this is required
project(example_project)            # this sets the project name

set( GSL_FOUND OFF )
set( GSL_CBLAS_FOUND OFF )

# Windows, but not for Cygwin and MSys where gsl-config is available
if( WIN32 AND NOT CYGWIN AND NOT MSYS )
  # look for headers
  find_path( GSL_INCLUDE_DIR
    NAMES gsl/gsl_cdf.h gsl/gsl_randist.h
    )
  if( GSL_INCLUDE_DIR )
    # look for gsl library
    find_library( GSL_LIBRARY
      NAMES gsl
    )
    if( GSL_LIBRARY )
      set( GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR} )
      get_filename_component( GSL_LIBRARY_DIRS ${GSL_LIBRARY} PATH )
      set( GSL_FOUND ON )
    endif( GSL_LIBRARY )

    # look for gsl cblas library
    find_library( GSL_CBLAS_LIBRARY
        NAMES gslcblas
      )
    if( GSL_CBLAS_LIBRARY )
      set( GSL_CBLAS_FOUND ON )
    endif( GSL_CBLAS_LIBRARY )

    set( GSL_LIBRARIES ${GSL_LIBRARY} ${GSL_CBLAS_LIBRARY} )
  endif( GSL_INCLUDE_DIR )

  mark_as_advanced(
    GSL_INCLUDE_DIR
    GSL_LIBRARY
    GSL_CBLAS_LIBRARY
  )
else( WIN32 AND NOT CYGWIN AND NOT MSYS )
  if( UNIX OR MSYS )
    find_program( GSL_CONFIG_EXECUTABLE gsl-config
      /usr/bin/
      /usr/local/bin
    )

    if( GSL_CONFIG_EXECUTABLE )
      set( GSL_FOUND ON )

      # run the gsl-config program to get cxxflags
      execute_process(
        COMMAND sh "${GSL_CONFIG_EXECUTABLE}" --cflags
        OUTPUT_VARIABLE GSL_CFLAGS
        RESULT_VARIABLE RET
        ERROR_QUIET
        )
      if( RET EQUAL 0 )
        string( STRIP "${GSL_CFLAGS}" GSL_CFLAGS )
        separate_arguments( GSL_CFLAGS )

        # parse definitions from cflags; drop -D* from CFLAGS
        string( REGEX MATCHALL "-D[^;]+"
          GSL_DEFINITIONS  "${GSL_CFLAGS}" )
        string( REGEX REPLACE "-D[^;]+;" ""
          GSL_CFLAGS "${GSL_CFLAGS}" )

        # parse include dirs from cflags; drop -I prefix
        string( REGEX MATCHALL "-I[^;]+"
          GSL_INCLUDE_DIRS "${GSL_CFLAGS}" )
        string( REPLACE "-I" ""
          GSL_INCLUDE_DIRS "${GSL_INCLUDE_DIRS}")
        string( REGEX REPLACE "-I[^;]+;" ""
          GSL_CFLAGS "${GSL_CFLAGS}")

        message("GSL_DEFINITIONS=${GSL_DEFINITIONS}")
        message("GSL_INCLUDE_DIRS=${GSL_INCLUDE_DIRS}")
        message("GSL_CFLAGS=${GSL_CFLAGS}")
      else( RET EQUAL 0 )
        set( GSL_FOUND FALSE )
      endif( RET EQUAL 0 )

      # run the gsl-config program to get the libs
      execute_process(
        COMMAND sh "${GSL_CONFIG_EXECUTABLE}" --libs
        OUTPUT_VARIABLE GSL_LIBRARIES
        RESULT_VARIABLE RET
        ERROR_QUIET
        )
      if( RET EQUAL 0 )
        string(STRIP "${GSL_LIBRARIES}" GSL_LIBRARIES )
        separate_arguments( GSL_LIBRARIES )

        # extract linkdirs (-L) for rpath (i.e., LINK_DIRECTORIES)
        string( REGEX MATCHALL "-L[^;]+"
          GSL_LIBRARY_DIRS "${GSL_LIBRARIES}" )
        string( REPLACE "-L" ""
          GSL_LIBRARY_DIRS "${GSL_LIBRARY_DIRS}" )
      else( RET EQUAL 0 )
        set( GSL_FOUND FALSE )
      endif( RET EQUAL 0 )

      MARK_AS_ADVANCED(
        GSL_CFLAGS
      )
      message( STATUS "Using GSL from ${GSL_PREFIX}" )
    else( GSL_CONFIG_EXECUTABLE )
      message( STATUS "FindGSL: gsl-config not found.")
    endif( GSL_CONFIG_EXECUTABLE )
  endif( UNIX OR MSYS )
endif( WIN32 AND NOT CYGWIN AND NOT MSYS )

if( GSL_FOUND )
  if( NOT GSL_FIND_QUIETLY )
    message( STATUS "FindGSL: Found both GSL headers and library" )
  endif( NOT GSL_FIND_QUIETLY )
else( GSL_FOUND )
  if( GSL_FIND_REQUIRED )
    message( FATAL_ERROR "FindGSL: Could not find GSL headers or library" )
  endif( GSL_FIND_REQUIRED )
endif( GSL_FOUND )



set(sources main.cpp cell.cpp rand.cpp types.cpp)

# You can use set(sources src/main.cpp) etc if you don't want to
# use globbing to find files automatically.

###############################################################################
## target definitions #########################################################
###############################################################################

# The data is just added to the executable, because in some IDEs (QtCreator)
# files are invisible when they are not explicitly part of the project.

# Just for example add some compiler flags.
set ( GCC_COVERAGE_COMPILE_FLAGS "-g -O -Wall  -pedantic  \
-pedantic-errors -Wextra -Waggregate-return -Wcast-align \
-Wcast-qual -Wconversion \
-Wdisabled-optimization \
-Werror -Wfloat-equal -Wformat=2 \
-Wformat-nonliteral -Wformat-security  \
-Wformat-y2k \
-Wimport  -Winit-self  -Winline \
-Winvalid-pch   \
-Wlong-long \
-Wmissing-field-initializers -Wmissing-format-attribute   \
-Wmissing-include-dirs -Wmissing-noreturn \
-Wpacked -Wpointer-arith \
-Wredundant-decls \
-Wshadow -Wstack-protector \
-Wstrict-aliasing=2 -Wswitch-default \
-Wswitch-enum \
-Wunreachable-code -Wunused \
-Wunused-parameter \
-Wvariadic-macros \
-Wwrite-strings")

#-Wpadded -Weffc++

set ( GCC_COVERAGE_COMPILE_FLAGS2 "-g   "  )
set ( GCC_COVERAGE_COMPILE_FLAGS3 "-g -O3 "  )
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS2})
set( CMAKE_VERBOSE_MAKEFILE on )

add_executable(testex ${sources})
target_link_libraries( testex ${GSL_LIBRARIES} )


