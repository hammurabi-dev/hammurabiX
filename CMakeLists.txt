cmake_minimum_required(VERSION 3.2.3)

#-------------- customized zone --------------#
set(CMAKE_CXX_COMPILER "g++")
option(ENABLE_TESTING "Enable testing for this project" OFF)
option(ENABLE_TIMING "Enable timing for this porect" ON)
option(ON_DOCKER "Build on docker image" ON)
option(BUILD_SHARED_LIBS "Build shared library" ON)
#-------------- instruction ------------------#
#
# make sure your compiler is GNU/INTEL, no CLANG, no way
# ENABLE_TESTING by default OFF, unless you are as paranoid as our developers
# ON_DOCKER by defauly ON, we highly recommend docker image for non-HPC tasks
# BUILD_SHARED_LIB by default ON, unless your IT department sucks
# 
# you have to specify your local paths of external libraries just below here
# in some special cases you have to modify find_path/find_library functions
# 
# update your CMake as much as you can
#
# if you add new modules/derived classes beyond original code
# please manually add source file paths to set(SRC_FILES ...) function
#
# complain any inconvenience to J.Wang (jiwang@sissa.it)
#
# @ enthusiastic or paranoid developers:
# we use GoogleTest for tests
# GoogleTest package is assembled INTO testing modules manually
# you can either install GoogleTest and cp src dir into install path
# or just download GoogleTest and specify root dir to GTEST_HINTS
#--------------------------------------------#

# we set GSL_HINTS just in case FindGSL fails
if(ON_DOCKER)
	set(GSL_HINTS /usr) 
	set(FFTW_HINTS /usr/local)
	set(CFITSIO_HINTS /usr/local)
	set(HEALPIX_HINTS /home/lab/Healpix_3.31/src/cxx/optimized_gcc)
	set(GTEST_HINTS /usr/local)
	set(CMAKE_INSTALL_PREFIX /usr/local/hammurabi)
else()
	message("apply customized fftw3/cfitsio/healpix paths")
	set(GSL_HINTS /u/shared/programs/x86_64/gsl/2.2/gnu/4.9.2)
	set(FFTW_HINTS /u/shared/programs/x86_64/fftw/3.3.4/gnu/4.9.2)
	set(CFITSIO_HINTS /home/jiwang/package/cfitsio)
	set(HEALPIX_HINTS /home/jiwang/package/Healpix_3.30/src/cxx/generic_gcc)
	set(GTEST_HINTS /home/jiwang/package/gtest)
	set(CMAKE_INSTALL_PREFIX /home/jiwang/package/hammurabi)
endif()
#---------------------------------------------#

project(hammurabiX CXX)
if(ENABLE_TESTING)
	enable_testing()
endif()

# setup cxx standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# compile flags under defferent options
if(ENABLE_TESTING AND ENABLE_TIMING)
	set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated")
elseif(ENABLE_TIMING)
	set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -DNDEBUG")
elseif(ENABLE_TESTING)
	set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -DNTIMING")
else()
	set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -DNTIMING -DNDEBUG")
endif()

# openmp and thread support
# if FindOpenMP fails, try add -fopenmp to CMAKE_CXX_FLAGS above
# the same solution applies to -pthread
include(FindOpenMP)
if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
	message(FATAL_ERROR "openmp unsupported")
endif()
if(ENABLE_TESTING)
	find_package(Threads REQUIRED)
	list(APPEND ALL_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif()

# we assemble include and external libs together
set(ALL_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include)
set(ALL_LIBRARIES)

# find sources
set(SRC_FILES 
	src/fields/gmf/brnd_jaffe.cc
	src/fields/gmf/breg_jaffe.cc
	src/fields/gmf/breg_verify.cc
	src/fields/gmf/breg_wmap.cc
	src/fields/gmf/breg.cc
	src/fields/gmf/brnd_es.cc
	src/fields/gmf/brnd_mhd.cc
	src/fields/gmf/brnd.cc
	src/fields/cre/cre.cc
	src/fields/cre/cre_num.cc
	src/fields/cre/cre_verify.cc
	src/fields/cre/cre_ana.cc
	src/fields/fe/fereg_verify.cc
	src/fields/fe/fereg.cc
	src/fields/fe/fernd_dft.cc
	src/fields/fe/fernd.cc
	src/fields/fe/fereg_ymw16.cc
	src/params/param.cc
	src/integrators/integrator.cc
	src/grids/grid.cc
	src/grids/grid_int.cc
	src/grids/grid_fernd.cc
	src/grids/grid_brnd.cc
	src/grids/grid_cre.cc
	src/grids/grid_fereg.cc
	src/grids/grid_breg.cc
	src/tools/tinyxml2.cc
	src/tools/namespace_toolkit.cc
)

# find HEALPix
find_path(HEALPIX_INCLUDE_DIR
	NAMES healpix_base.h
	HINTS 
	${HEALPIX_HINTS}/include
)
find_library(HEALPIX_CXX_LIBRARY
	NAMES healpix_cxx
	HINTS 
	${HEALPIX_HINTS}/lib
)
find_library(HEALPIX_SUP_LIBRARY
	NAMES cxxsupport
	HINTS 
	${HEALPIX_HINTS}/lib
)
find_library(HEALPIX_UTI_LIBRARY
	NAMES c_utils
	HINTS 
	${HEALPIX_HINTS}/lib
)
find_library(HEALPIX_FFT_LIBRARY
	NAMES fftpack
	HINTS 
	${HEALPIX_HINTS}/lib
)
find_library(HEALPIX_SHA_LIBRARY
	NAMES sharp
	HINTS
	${HEALPIX_HINTS}/lib
)
if(HEALPIX_INCLUDE_DIR)
	list(APPEND ALL_INCLUDE_DIR ${HEALPIX_INCLUDE_DIR})
else()
	message(FATAL_ERROR "${HEALPIX_INCLUDE_DIR}")
endif()
if(HEALPIX_CXX_LIBRARY)
	list(APPEND ALL_LIBRARIES ${HEALPIX_CXX_LIBRARY})
else()
	message(FATAL_ERROR "${HEALPIX_CXX_LIBRARY}")
endif()
if(HEALPIX_SUP_LIBRARY)
	list(APPEND ALL_LIBRARIES ${HEALPIX_SUP_LIBRARY})
else()
	message(FATAL_ERROR "${HEALPIX_SUP_LIBRARY}")
endif()
if(HEALPIX_UTI_LIBRARY)
	list(APPEND ALL_LIBRARIES ${HEALPIX_UTI_LIBRARY})
else()
	message(FATAL_ERROR "${HEALPIX_UTI_LIBRARY}")
endif()
if(HEALPIX_FFT_LIBRARY)
	list(APPEND ALL_LIBRARIES ${HEALPIX_FFT_LIBRARY})
else()
	message(FATAL_ERROR "${HEALPIX_FFT_LIBRARY}")
endif()
if(HEALPIX_SHA_LIBRARY)
	list(APPEND ALL_LIBRARIES ${HEALPIX_SHA_LIBRARY})
else()
	message(FATAL_ERROR "${HEALPIX_SHA_LIBRARY}")
endif()


# find CFITSIO
find_path(CFITSIO_INCLUDE_DIR
	NAMES fitsio.h
	HINTS
	${CFITSIO_HINTS}/include
)
find_library(CFITSIO_LIBRARY
	NAMES cfitsio
	HINTS
	${CFITSIO_HINTS}/lib
)
if(CFITSIO_INCLUDE_DIR)
	list(APPEND ALL_INCLUDE_DIR ${CFITSIO_INCLUDE_DIR})
else()
	message(FATAL_ERROR "${CFITSIO_INCLUDE_DIR}")
endif()
if(CFITSIO_LIBRARY)
	list(APPEND ALL_LIBRARIES ${CFITSIO_LIBRARY})
else()
	message(FATAL_ERROR "${CFITSIO_LIBRARY}")
endif()

# find FFTW, FFTW_OMP
find_path(FFTW_INCLUDE_DIR
	NAMES fftw3.h
	HINTS
	${FFTW_HINTS}/include
)
find_library(FFTW_LIBRARY
	NAMES fftw3
	HINTS
	${FFTW_HINTS}/lib	
)
find_library(FFTW_OMP_LIBRARY
	NAMES fftw3_omp
	HINTS
	${FFTW_HINTS}/lib
)
if(FFTW_INCLUDE_DIR)
	list(APPEND ALL_INCLUDE_DIR ${FFTW_INCLUDE_DIR})
else()
	message(FATAL_ERROR "${FFTW_INCLUDE_DIR}")
endif()
if(FFTW_LIBRARY)
	list(APPEND ALL_LIBRARIES ${FFTW_LIBRARY})
else()
	message(FATAL_ERROR "${FFTW_LIBRARY}")
endif()
if(FFTW_OMP_LIBRARY)
	list(APPEND ALL_LIBRARIES ${FFTW_OMP_LIBRARY})
else()
	message(FATAL_ERROR "${FFTW_OMP_LIBRARY}")
endif()

# find GSL, GSL-cblas
# if FindGSL fails, try find_path & find_library functions manually
find_package(GSL REQUIRED)

# find GSL, GSL-cblas manually
if(NOT GSL_FOUND)
find_path(GSL_INCLUDE_DIR
	NAMES gsl_math.h
	HINTS 
	${GSL_HINTS}/include
	PATH_SUFFIXES gsl
)
find_library(GSL_LIBRARY
	NAMES gsl
	HINTS 
	${GSL_HINTS}/lib
	PATH_SUFFIXES x86_64-linux-gnu
)
find_library(GSL_CBLAS_LIBRARY
	NAMES gslcblas
	HINTS
	${GSL_HINTS}/lib
	PATH_SUFFIXES x86_64-linux-gnu
)
if(GSL_INCLUDE_DIR)
	list(APPEND ALL_INCLUDE_DIR ${GSL_INCLUDE_DIR})
else()
	message(FATAL_ERROR "${GSL_INCLUDE_DIR}")
endif()
if(GSL_LIBRARY)
	list(APPEND ALL_LIBRARIES ${GSL_LIBRARY})
else()
	message(FATAL_ERROR "${GSL_LIBRARY}")
endif()
if(GSL_CBLAS_LIBRARY)
	list(APPEND ALL_LIBRARIES ${GSL_CBLAS_LIBRARY})
else()
	message(FATAL_ERROR "${GSL_CBLAS_LIBRARY}")
endif()
endif()

# find GoogleTest
# GoogleTest package must be assembled INTO testing modules manually
# you can either install GoogleTest and cp src dir into install path
# or just download GoogleTest
if(ENABLE_TESTING)
	find_path(GTEST_INC_DIR 
		NAMES gtest/gtest.h
  		HINTS
		${GTEST_HINTS}/include
		${GTEST_HINTS}/googletest/include
	)
	if(NOT GTEST_INC_DIR)
  		message(FATAL_ERROR ${GTEST_INC_DIR})
	endif()
	find_path(GTEST_SOURCE_DIR gtest-all.cc
  		${GTEST_HINTS}/src
  		${GTEST_HINTS}/googletest/src
	)
	if(NOT GTEST_SOURCE_DIR)
  		message(FATAL_ERROR ${GTEST_SOURCE_DIR})
	endif()
	set(GTEST_LIB_SOURCES ${GTEST_SOURCE_DIR}/gtest-all.cc)
	set(GTEST_MAIN_SOURCES ${GTEST_SOURCE_DIR}/gtest_main.cc)
	set(GTEST_INCLUDE_DIRS ${GTEST_INC_DIR} ${GTEST_SOURCE_DIR} ${GTEST_SOURCE_DIR}/..)
endif()

# assemble lib
include_directories(${ALL_INCLUDE_DIR})
add_library(hammurabi ${SRC_FILES})
target_link_libraries(hammurabi ${ALL_LIBRARIES} GSL::gsl GSL::gslcblas)

# build testing cases
if(ENABLE_TESTING)
	add_subdirectory(tests)
endif()

# build executable
add_executable(hamx src/main/hammurabi.cc)
target_link_libraries(hamx hammurabi)
# copy template parameter file into build directory
file(COPY ${CMAKE_CURRENT_LIST_DIR}/templates/params_template.xml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# install package
install(TARGETS hamx DESTINATION bin)
install(FILES 
	include/breg.h
	include/brnd.h
	include/cgs_units_file.h
	include/cre.h
	include/fereg.h
	include/fernd.h
	include/integrator.h
	include/tinyxml2.h
	include/timer.h
	include/namespace_toolkit.h
	include/param.h
	include/grid.h	
	DESTINATION include
)
install(TARGETS hammurabi DESTINATION lib)
message("hammurabi X will be installed at ${CMAKE_INSTALL_PREFIX}")
