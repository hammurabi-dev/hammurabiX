CMAKE_MINIMUM_REQUIRED(VERSION 3.2.3)

#-------------- customized zone --------------#
SET(CMAKE_CXX_COMPILER "g++")
OPTION(ENABLE_TESTING "Enable testing for this project" ON)
OPTION(ENABLE_TIMING "Enable timing for this porect" ON)
OPTION(ON_DOCKER "Build on docker image" ON)
OPTION(BUILD_SHARED_LIBS "Build shared library" ON)
OPTION(ENABLE_REPORT "Enable verbose report" ON)
#-------------- instruction ------------------#
#
# make sure your compiler is GNU or INTEL,
# ENABLE_TESTING by default OFF, unless you are as paranoid as us
# ENABLE_TIMING is also not necessary for users
# ON_DOCKER by defauly ON, we highly recommend Docker image for non-HPC tasks
# but if install manually, switch it off and modify the path hints listed below
# BUILD_SHARED_LIB by default ON
# you will be overwhelmed by ENABLE_REPORT, switch it off for non-testing tasks
# 
# you have to specify your local paths of external libraries just below here
# in some special cases you have to modify FIND_PATH/FIND_LIBRARY functions
#
# in some special cases, you may want to manually fix LFLAGS or CFLAGS
# 
# update your CMake as much as you can
#
# if you add new modules/derived classes beyond original code
# please manually add source file paths to SET(SRC_FILES ...) function
#
# the last resort would be manually calling homemade Makefile building system
# you can find cached building files in mfcache
#
# if even the last resort fails you
# complain your issue to J.Wang (jiwang@sissa.it)
#
# @ enthusiastic or paranoid developers:
# we use Google Test for tests
# Google Test package is assembled INTO testing modules manually
# you can either install GoogleTest and cp src dir into install path
# or just download GoogleTest and specify root dir to GTEST_HINTS
#--------------------------------------------#
# we SET GSL_HINTS just in case FindGSL fails
IF(ON_DOCKER)
	SET(GSL_HINTS /usr) 
	SET(FFTW_HINTS /usr/local)
	SET(CFITSIO_HINTS /usr/local)
	SET(HEALPIX_HINTS /home/lab/healpix/src/cxx/optimized_gcc)
	SET(GTEST_HINTS /usr/local)
	SET(INSTALL_ROOT_DIR /usr/local/hammurabi)
ELSE()
	MESSAGE("apply customized fftw3/cfitsio/healpix paths")
	SET(GSL_HINTS $HOME/local/gsl)
	SET(FFTW_HINTS $HOME/local/fftw)
	SET(CFITSIO_HINTS $HOME/local/cfitsio)
	SET(HEALPIX_HINTS $HOME/local/Healpix_3.50/src/cxx/optimized_gcc)
	SET(GTEST_HINTS $HOME/local/googletest)
	SET(INSTALL_ROOT_DIR $HOME/local/hammurabi)
	SET(CMAKE_EXE_LINKER_FLAGS -lcurl)
ENDIF()
#---------------------------------------------#

PROJECT(hammurabiX CXX)

# SETup cxx standard
SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_STANDARD_REQUIRED YES)

# compile flags under defferent options
IF(ENABLE_TESTING AND ENABLE_TIMING)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated")
ELSEIF(ENABLE_TIMING)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -DNDEBUG")
ELSEIF(ENABLE_TESTING)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -DNTIMING")
ELSE()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -DNTIMING -DNDEBUG")
ENDIF()

IF (ENABLE_REPORT)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE")
ENDIF()

# openmp and thread support
# if FindOpenMP fails, try add -fopenmp to CMAKE_CXX_FLAGS above
# the same solution applies to -pthread
include(FindOpenMP)
IF(OPENMP_FOUND)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ELSE()
	MESSAGE(FATAL_ERROR "openmp unsupported")
ENDIF()
IF(ENABLE_TESTING)
	FIND_PACKAGE(Threads REQUIRED)
	LIST(APPEND ALL_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
ENDIF()

# we assemble include and external libs together
SET(ALL_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
SET(ALL_LIBRARIES)

# find sources
SET(SRC_FILES
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/breg.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/breg_jaffe.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/breg_lsa.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/breg_unif.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/brnd.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/brnd_es.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/brnd_mhd.cc

	${CMAKE_CURRENT_LIST_DIR}/source/field/cre/cre.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/cre/cre_ana.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/cre/cre_unif.cc

	${CMAKE_CURRENT_LIST_DIR}/source/field/te/tereg.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/te/tereg_unif.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/te/tereg_ymw16.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/te/ternd.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/te/ternd_dft.cc

	${CMAKE_CURRENT_LIST_DIR}/source/param/param.cc

	${CMAKE_CURRENT_LIST_DIR}/source/integrator/integrator.cc

	${CMAKE_CURRENT_LIST_DIR}/source/grid/grid.cc
	${CMAKE_CURRENT_LIST_DIR}/source/grid/grid_obs.cc
	${CMAKE_CURRENT_LIST_DIR}/source/grid/grid_tereg.cc
	${CMAKE_CURRENT_LIST_DIR}/source/grid/grid_breg.cc
	${CMAKE_CURRENT_LIST_DIR}/source/grid/grid_ternd.cc
	${CMAKE_CURRENT_LIST_DIR}/source/grid/grid_brnd.cc
	${CMAKE_CURRENT_LIST_DIR}/source/grid/grid_cre.cc

	${CMAKE_CURRENT_LIST_DIR}/source/tools/tinyxml2.cc
	${CMAKE_CURRENT_LIST_DIR}/source/tools/toolkit.cc
)

# find HEALPix
FIND_PATH(HEALPIX_INCLUDE_DIR
	NAMES healpix_base.h
	HINTS 
	${HEALPIX_HINTS}/include
)
FIND_LIBRARY(HEALPIX_CXX_LIBRARY
	NAMES healpix_cxx
	HINTS 
	${HEALPIX_HINTS}/lib
)
FIND_LIBRARY(HEALPIX_SUP_LIBRARY
	NAMES cxxsupport
	HINTS 
	${HEALPIX_HINTS}/lib
)
FIND_LIBRARY(HEALPIX_UTI_LIBRARY
	NAMES c_utils
	HINTS 
	${HEALPIX_HINTS}/lib
)
FIND_LIBRARY(HEALPIX_FFT_LIBRARY
	NAMES pocketfft
	HINTS 
	${HEALPIX_HINTS}/lib
)
FIND_LIBRARY(HEALPIX_SHA_LIBRARY
	NAMES sharp
	HINTS
	${HEALPIX_HINTS}/lib
)
IF(HEALPIX_INCLUDE_DIR)
	LIST(APPEND ALL_INCLUDE_DIR ${HEALPIX_INCLUDE_DIR})
ELSE()
	MESSAGE(FATAL_ERROR "${HEALPIX_INCLUDE_DIR}")
ENDIF()
IF(HEALPIX_CXX_LIBRARY)
	LIST(APPEND ALL_LIBRARIES ${HEALPIX_CXX_LIBRARY})
ELSE()
	MESSAGE(FATAL_ERROR "${HEALPIX_CXX_LIBRARY}")
ENDIF()
IF(HEALPIX_SUP_LIBRARY)
	LIST(APPEND ALL_LIBRARIES ${HEALPIX_SUP_LIBRARY})
ELSE()
	MESSAGE(FATAL_ERROR "${HEALPIX_SUP_LIBRARY}")
ENDIF()
IF(HEALPIX_UTI_LIBRARY)
	LIST(APPEND ALL_LIBRARIES ${HEALPIX_UTI_LIBRARY})
ELSE()
	MESSAGE(FATAL_ERROR "${HEALPIX_UTI_LIBRARY}")
ENDIF()
IF(HEALPIX_FFT_LIBRARY)
	LIST(APPEND ALL_LIBRARIES ${HEALPIX_FFT_LIBRARY})
ELSE()
	MESSAGE(FATAL_ERROR "${HEALPIX_FFT_LIBRARY}")
ENDIF()
IF(HEALPIX_SHA_LIBRARY)
	LIST(APPEND ALL_LIBRARIES ${HEALPIX_SHA_LIBRARY})
ELSE()
	MESSAGE(FATAL_ERROR "${HEALPIX_SHA_LIBRARY}")
ENDIF()


# find CFITSIO
FIND_PATH(CFITSIO_INCLUDE_DIR
	NAMES fitsio.h
	HINTS
	${CFITSIO_HINTS}/include
)
FIND_LIBRARY(CFITSIO_LIBRARY
	NAMES cfitsio
	HINTS
	${CFITSIO_HINTS}/lib
)
IF(CFITSIO_INCLUDE_DIR)
	LIST(APPEND ALL_INCLUDE_DIR ${CFITSIO_INCLUDE_DIR})
ELSE()
	MESSAGE(FATAL_ERROR "${CFITSIO_INCLUDE_DIR}")
ENDIF()
IF(CFITSIO_LIBRARY)
	LIST(APPEND ALL_LIBRARIES ${CFITSIO_LIBRARY})
ELSE()
	MESSAGE(FATAL_ERROR "${CFITSIO_LIBRARY}")
ENDIF()

# find FFTW, FFTW_OMP
FIND_PATH(FFTW_INCLUDE_DIR
	NAMES fftw3.h
	HINTS
	${FFTW_HINTS}/include
)
FIND_LIBRARY(FFTW_LIBRARY
	NAMES fftw3
	HINTS
	${FFTW_HINTS}/lib	
)
FIND_LIBRARY(FFTW_OMP_LIBRARY
	NAMES fftw3_omp
	HINTS
	${FFTW_HINTS}/lib
)
IF(FFTW_INCLUDE_DIR)
	LIST(APPEND ALL_INCLUDE_DIR ${FFTW_INCLUDE_DIR})
ELSE()
	MESSAGE(FATAL_ERROR "${FFTW_INCLUDE_DIR}")
ENDIF()
IF(FFTW_LIBRARY)
	LIST(APPEND ALL_LIBRARIES ${FFTW_LIBRARY})
ELSE()
	MESSAGE(FATAL_ERROR "${FFTW_LIBRARY}")
ENDIF()
IF(FFTW_OMP_LIBRARY)
	LIST(APPEND ALL_LIBRARIES ${FFTW_OMP_LIBRARY})
ELSE()
	MESSAGE(FATAL_ERROR "${FFTW_OMP_LIBRARY}")
ENDIF()

# find GSL, GSL-cblas
# IF FindGSL fails, try FIND_PATH & FIND_LIBRARY functions manually
FIND_PACKAGE(GSL REQUIRED)

# find GSL, GSL-cblas manually
IF(NOT GSL_FOUND)
FIND_PATH(GSL_INCLUDE_DIR
	NAMES gsl_math.h
	HINTS 
	${GSL_HINTS}/include
	PATH_SUFFIXES gsl
)
FIND_LIBRARY(GSL_LIBRARY
	NAMES gsl
	HINTS 
	${GSL_HINTS}/lib
	PATH_SUFFIXES x86_64-linux-gnu
)
FIND_LIBRARY(GSL_CBLAS_LIBRARY
	NAMES gslcblas
	HINTS
	${GSL_HINTS}/lib
	PATH_SUFFIXES x86_64-linux-gnu
)
IF(GSL_INCLUDE_DIR)
	LIST(APPEND ALL_INCLUDE_DIR ${GSL_INCLUDE_DIR})
ELSE()
	MESSAGE(FATAL_ERROR "${GSL_INCLUDE_DIR}")
ENDIF()
IF(GSL_LIBRARY)
	LIST(APPEND ALL_LIBRARIES ${GSL_LIBRARY})
ELSE()
	MESSAGE(FATAL_ERROR "${GSL_LIBRARY}")
ENDIF()
IF(GSL_CBLAS_LIBRARY)
	LIST(APPEND ALL_LIBRARIES ${GSL_CBLAS_LIBRARY})
ELSE()
	MESSAGE(FATAL_ERROR "${GSL_CBLAS_LIBRARY}")
ENDIF()
ENDIF()

# find GoogleTest
# GoogleTest package must be assembled INTO testing modules manually
# you can either install GoogleTest and cp src dir into install path
# or just download GoogleTest
IF(ENABLE_TESTING)
	FIND_PATH(GTEST_INC_DIR 
		NAMES gtest/gtest.h
  		HINTS
		${GTEST_HINTS}/include
		${GTEST_HINTS}/googletest/include
	)
	IF(NOT GTEST_INC_DIR)
  		MESSAGE(FATAL_ERROR ${GTEST_INC_DIR})
	ENDIF()
	FIND_PATH(GTEST_SOURCE_DIR gtest-all.cc
  		${GTEST_HINTS}/src
  		${GTEST_HINTS}/googletest/src
	)
	IF(NOT GTEST_SOURCE_DIR)
  		MESSAGE(FATAL_ERROR ${GTEST_SOURCE_DIR})
	ENDIF()
	SET(GTEST_LIB_SOURCES ${GTEST_SOURCE_DIR}/gtest-all.cc)
	SET(GTEST_MAIN_SOURCES ${GTEST_SOURCE_DIR}/gtest_main.cc)
	SET(GTEST_INCLUDE_DIRS ${GTEST_INC_DIR} ${GTEST_SOURCE_DIR} ${GTEST_SOURCE_DIR}/..)
ENDIF()

# assemble lib
INCLUDE_DIRECTORIES(${ALL_INCLUDE_DIR})
ADD_LIBRARY(hammurabi ${SRC_FILES})
TARGET_LINK_LIBRARIES(hammurabi ${ALL_LIBRARIES} GSL::gsl GSL::gslcblas)

# build testing cases
IF(ENABLE_TESTING)
	ENABLE_TESTING()
	ADD_SUBDIRECTORY(tests)
ENDIF()

# build executable
ADD_EXECUTABLE(hamx source/main/hammurabi.cc)
TARGET_LINK_LIBRARIES(hamx hammurabi)
# copy template parameter file into build directory
FILE(COPY
	${CMAKE_CURRENT_LIST_DIR}/templates/params_template.xml
DESTINATION
	${CMAKE_CURRENT_BINARY_DIR}
)

# install package
SET(CMAKE_INSTALL_PREFIX ${INSTALL_ROOT_DIR})
INSTALL(TARGETS hamx DESTINATION bin)
INSTALL(FILES
	${CMAKE_CURRENT_LIST_DIR}/include/cgs_units.h
	${CMAKE_CURRENT_LIST_DIR}/include/bfield.h
	${CMAKE_CURRENT_LIST_DIR}/include/crefield.h
	${CMAKE_CURRENT_LIST_DIR}/include/tefield.h
	${CMAKE_CURRENT_LIST_DIR}/include/integrator.h
	${CMAKE_CURRENT_LIST_DIR}/include/tinyxml2.h
	${CMAKE_CURRENT_LIST_DIR}/include/timer.h
	${CMAKE_CURRENT_LIST_DIR}/include/hamvec.h
	${CMAKE_CURRENT_LIST_DIR}/include/hampixp.h
	${CMAKE_CURRENT_LIST_DIR}/include/toolkit.h
	${CMAKE_CURRENT_LIST_DIR}/include/param.h
	${CMAKE_CURRENT_LIST_DIR}/include/grid.h	
	DESTINATION include
)
INSTALL(TARGETS hammurabi DESTINATION lib)
MESSAGE("hammurabi X will be installed at ${CMAKE_INSTALL_PREFIX}")
