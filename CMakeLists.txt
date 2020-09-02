CMAKE_MINIMUM_REQUIRED(VERSION 3.2.3)

#-------------- customized zone --------------#

SET(CMAKE_CXX_COMPILER "g++")
OPTION(ENABLE_TESTING "Enable testing for this project" ON)
OPTION(ENABLE_TIMING "Enable timing for this porect" ON)
OPTION(ON_DOCKER "Build on docker image" ON)
OPTION(BUILD_SHARED_LIBS "Build shared library" ON)
OPTION(ENABLE_REPORT "Enable verbose report" ON)

#-------------- instruction ------------------#

# ENABLE_TESTING by default ON, for package build testing,
# ENABLE_TIMING is for checking performance,
# ON_DOCKER by defauly ON, 
# we highly recommend Docker image for non-HPC tasks,
# but if install manually,
# switch it off and modify the path hints LISTED BELOW,
# BUILD_SHARED_LIB by default ON,
# you will be overwhelmed by ENABLE_REPORT,
# switch it off for non-testing tasks,
# 
# you have to specify your local paths of external libraries just below here,
# in some special cases you have to modify FIND_PATH/FIND_LIBRARY functions,
#
# in some special cases, you may want to manually fix LFLAGS or CFLAGS,
#
# if you add new modules/derived classes beyond original code,
# please manually add source file paths to SET(SRC_FILES ...) function,
#
# the last resort would be manually calling homemade Makefile building system,
# you can find cached building files in "cache",
#
# we use Google Test for tests,
# Google Test package is assembled INTO testing modules manually,
# you can either install GoogleTest and cp src dir into install path,
# or just download GoogleTest and specify root dir to GTEST_HINTS

#--------------------------------------------#

# we SET GSL_HINTS just in case FindGSL fails
IF(ON_DOCKER)
	SET(GSL_HINTS /usr) 
	SET(FFTW_HINTS /tmp/local)
	SET(GTEST_HINTS /tmp/local)
	SET(INSTALL_ROOT_DIR /tmp/local/hammurabi)
ELSE()
	MESSAGE("apply GSL FFTW3 GoogleTest paths")
	SET(GSL_HINTS $HOME/local/gsl)
	SET(FFTW_HINTS $HOME/local/fftw)
	SET(GTEST_HINTS $HOME/local/gtest)
	SET(INSTALL_ROOT_DIR $HOME/local/hammurabi)
	SET(CMAKE_EXE_LINKER_FLAGS)
ENDIF()

#---------------------------------------------#

PROJECT(hammurabiX CXX)

# setup cxx standard

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_STANDARD_REQUIRED YES)

# compile flags under defferent options

IF(ENABLE_TESTING AND ENABLE_TIMING)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -g")
ELSEIF(ENABLE_TIMING)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -g -DNDEBUG")
ELSEIF(ENABLE_TESTING)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -g -DNTIMING")
ELSE()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -g -DNTIMING -DNDEBUG")
ENDIF()

IF (ENABLE_REPORT)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE")
ENDIF()

# openmp and thread support
# if FindOpenMP fails, try add -fopenmp to CMAKE_CXX_FLAGS above
# the same solution applies to -pthread

INCLUDE(FindOpenMP)
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
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/breg_cart.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/breg_helix.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/breg_jf12.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/brnd.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/brnd_es.cc
	${CMAKE_CURRENT_LIST_DIR}/source/field/b/brnd_jf12.cc
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

	${CMAKE_CURRENT_LIST_DIR}/source/xml/tinyxml2.cc

	${CMAKE_CURRENT_LIST_DIR}/source/pipeline/pipeline.cc
)

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

# assemble library

INCLUDE_DIRECTORIES(${ALL_INCLUDE_DIR})
ADD_LIBRARY(hammurabi ${SRC_FILES})
TARGET_LINK_LIBRARIES(hammurabi ${ALL_LIBRARIES} GSL::gsl GSL::gslcblas)

# build testing cases

IF(ENABLE_TESTING)
	ENABLE_TESTING()
	ADD_SUBDIRECTORY(tests)
ENDIF()

# build executable

ADD_EXECUTABLE(hamx source/main/main_std.cc)
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
	${CMAKE_CURRENT_LIST_DIR}/include/tinyxml2.h
	${CMAKE_CURRENT_LIST_DIR}/include/hamtype.h
	${CMAKE_CURRENT_LIST_DIR}/include/hamunits.h
	${CMAKE_CURRENT_LIST_DIR}/include/hamvec.h
	${CMAKE_CURRENT_LIST_DIR}/include/hamp.h
	${CMAKE_CURRENT_LIST_DIR}/include/hamdis.h
	${CMAKE_CURRENT_LIST_DIR}/include/hamio.h
	${CMAKE_CURRENT_LIST_DIR}/include/hamsk.h
	${CMAKE_CURRENT_LIST_DIR}/include/toolkit.h
	${CMAKE_CURRENT_LIST_DIR}/include/timer.h
	${CMAKE_CURRENT_LIST_DIR}/include/bfield.h
	${CMAKE_CURRENT_LIST_DIR}/include/crefield.h
	${CMAKE_CURRENT_LIST_DIR}/include/tefield.h
	${CMAKE_CURRENT_LIST_DIR}/include/integrator.h
	${CMAKE_CURRENT_LIST_DIR}/include/param.h
	${CMAKE_CURRENT_LIST_DIR}/include/grid.h	
	DESTINATION include
)
INSTALL(TARGETS hammurabi DESTINATION lib)
MESSAGE("hammurabi X will be installed at ${CMAKE_INSTALL_PREFIX}")
