cmake_minimum_required(VERSION 2.8)
project(hammurabiX CXX)

#-------------- customized zone --------------#
option(ENABLE_TESTING "Enable testing for this project" OFF)
option(ENABLE_TIMING "Enable timing for this porect" ON)
option(ON_DOCKER "Build on docker image" OFF)
option(BUILD_SHARED_LIBS "Build shared library" OFF)
#-------------- instruction ------------------#
# you have to specify your local paths of external libraries just below here 
# in some special cases you have to modify find_path/find_library functions
#
# if you add new modules/derived classes beyond original code
# please manually add source file paths to set(SRC_FILES ...) function
#

if(ON_DOCKER)
	set(GSL_HINTS /usr)
	set(FFTW_HINTS /usr/local)
	set(CFITSIO_HINTS /usr/local)
	set(HEALPIX_HINTS /home/lab/Healpix_3.31/src/cxx/optimized_gcc)
	set(CMAKE_INSTALL_PREFIX /usr/local/hammurabi)
else()
	message("apply customized gsl/fftw3/cfitsio/healpix paths")
	set(GSL_HINTS /u/shared/programs/x86_64/gsl/2.2/gnu/4.9.2)
	set(FFTW_HINTS /u/shared/programs/x86_64/fftw/3.3.4/gnu/4.9.2)
	set(CFITSIO_HINTS /home/jiwang/package/cfitsio)
	set(HEALPIX_HINTS /home/jiwang/package/Healpix_3.30/src/cxx/generic_gcc)
	set(CMAKE_INSTALL_PREFIX /home/jiwang/package/hammurabi)
endif()

#---------------------------------------------#

# setup cxx standard
set(CMAKE_CXX_STANDARD 14)

# compile flags under defferent options
if(ENABLE_TESTING AND ENABLE_TIMING)
	set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -fopenmp -pthread")
elseif(ENABLE_TIMING)
	set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -fopenmp -pthread -DNDEBUG")
elseif(ENABLE_TESTING)
	set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -fopenmp -pthread -DNTIMING")
else()
	set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -fopenmp -pthread -DNTIMING -DNDEBUG")
endif()

set(ALL_INCLUDE_DIR include)
set(ALL_LIBRARIES)

# find sources
set(SRC_FILES 
	src/fields/gmf/brnd_jaffe.cc
	src/fields/gmf/breg_jaffe.cc
	src/fields/gmf/breg_verify.cc
	src/fields/gmf/breg_wmap.cc
	src/fields/gmf/breg.cc
	src/fields/gmf/brnd_global.cc
	src/fields/gmf/brnd_local.cc
	src/fields/gmf/brnd.cc
	src/fields/cre/cre.cc
	src/fields/cre/cre_num.cc
	src/fields/cre/cre_verify.cc
	src/fields/cre/cre_ana.cc
	src/fields/fe/fereg_verify.cc
	src/fields/fe/fereg.cc
	src/fields/fe/fernd_global.cc
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
	message(FATAL_ERROR "${HEALPIX_INCLUDE_DIR} not found")
endif()
if(HEALPIX_CXX_LIBRARY)
	list(APPEND ALL_LIBRARIES ${HEALPIX_CXX_LIBRARY})
else()
	message(FATAL_ERROR "${HEALPIX_CXX_LIBRARY} not found")
endif()
if(HEALPIX_SUP_LIBRARY)
	list(APPEND ALL_LIBRARIES ${HEALPIX_SUP_LIBRARY})
else()
	message(FATAL_ERROR "${HEALPIX_SUP_LIBRARY} not found")
endif()
if(HEALPIX_UTI_LIBRARY)
	list(APPEND ALL_LIBRARIES ${HEALPIX_UTI_LIBRARY})
else()
	message(FATAL_ERROR "${HEALPIX_UTI_LIBRARY} not found")
endif()
if(HEALPIX_FFT_LIBRARY)
	list(APPEND ALL_LIBRARIES ${HEALPIX_FFT_LIBRARY})
else()
	message(FATAL_ERROR "${HEALPIX_FFT_LIBRARY} not found")
endif()
if(HEALPIX_SHA_LIBRARY)
	list(APPEND ALL_LIBRARIES ${HEALPIX_SHA_LIBRARY})
else()
	message(FATAL_ERROR "${HEALPIX_SHA_LIBRARY} not found")
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
	message(FATAL_ERROR "${CFITSIO_INCLUDE_DIR} not found")
endif()
if(CFITSIO_LIBRARY)
	list(APPEND ALL_LIBRARIES ${CFITSIO_LIBRARY})
else()
	message(FATAL_ERROR "${CFITSIO_LIBRARY} not found")
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
	message(FATAL_ERROR "${FFTW_INCLUDE_DIR} not found")
endif()
if(FFTW_LIBRARY)
	list(APPEND ALL_LIBRARIES ${FFTW_LIBRARY})
else()
	message(FATAL_ERROR "${FFTW_LIBRARY} not found")
endif()
if(FFTW_OMP_LIBRARY)
	list(APPEND ALL_LIBRARIES ${FFTW_OMP_LIBRARY})
else()
	message(FATAL_ERROR "${FFTW_OMP_LIBRARY} not found")
endif()

# find GSL, GSL_CBLAS
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
	message(FATAL_ERROR "${GSL_INCLUDE_DIR} not found")
endif()
if(GSL_LIBRARY)
	list(APPEND ALL_LIBRARIES ${GSL_LIBRARY})
else()
	message(FATAL_ERROR "${GSL_LIBRARY} not found")
endif()
if(GSL_CBLAS_LIBRARY)
	list(APPEND ALL_LIBRARIES ${GSL_CBLAS_LIBRARY})
else()
	message(FATAL_ERROR "${GSL_CBLAS_LIBRARY} not found")
endif()

# debug info
#message("suffix: ${CMAKE_FIND_LIBRARY_SUFFIXES}")

# assemble static lib
include_directories(${ALL_INCLUDE_DIR})
add_library(hammurabi ${SRC_FILES})
target_link_libraries(hammurabi ${ALL_LIBRARIES})

# build executable
add_executable(hamx src/main/hammurabi.cc)
target_link_libraries(hamx hammurabi)

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