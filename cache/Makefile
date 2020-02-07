# prepared just in case CMake meets any issue on a HPC server

ifeq ($(CXX),)
CXX=g++
endif

ifeq ($(GSL),)
GSL=/usr/local/gsl
endif

ifeq ($(FFTW),)
FFTW=/usr/local/fftw
endif

ifeq ($(CFITSIO),)
CFITSIO=/usr/local/cifitsio
endif

ifeq ($(HEALPIX),)
HEALPIX=/usr/local/healpix
endif

ifeq ($(EXE),)
EXE=hamx
endif

SRCDIR=../source
OBJDIR=../build
SOURCES=$(shell find $(SRCDIR) -type f -name "*.cc")
OBJ=$(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(SOURCES:.cc=.o))
INCDIR=../include
HEADS=$(wildcard $(INCDIR)*.h)
LIBDIR=../lib
EXEDIR=../bin

# COMPILE OPTIONS
CXXFLAGS = -std=c++14 -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -fopenmp -pthread $(DEBUG)

EXT_INC = -I$(CFITSIO)/include -I$(GSL)/include -I$(FFTW)/include \
	-I$(HEALPIX)/include/healpix_cxx \
	-I $(INCDIR) -I.
CXXFLAGS += $(EXT_INC)
 
ALL_L = -L$(CFITSIO)/lib -L$(GSL)/lib -L$(FFTW)/lib \
	-L$(HEALPIX)/lib
LD_FLAGS = $(ALL_L) \
	-lhealpix_cxx -lcfitsio \
	-lgsl -lgslcblas -lm -lfftw3 -lfftw3_omp -lfftw3_threads -fopenmp $(DEBUG)
LFLAGS = $(LD_FLAGS)

# RULES

default: $(EXE)

$(EXE): | $(EXEDIR)

$(EXEDIR):
	mkdir -p $(EXEDIR)

$(EXE): $(OBJ)
	$(CXX) $^ -o $(EXEDIR)/$(EXE) $(LFLAGS)
	cp ../templates/* $(EXEDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(HEADS)
	mkdir -p $(dir $@)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(LIBDIR)/libhammurabi.so: $(OBJ)
	mkdir -p $(LIBDIR)
	$(CXX) -o $@ -shared -fPIC $(LD_FLAGS) $^

help:
	echo $(OBJ)

lib: $(LIBDIR)/libhammurabi.so

.PHONY: clean
clean:
	rm -rf $(OBJDIR)

.PHONY: wipe
wipe:
	rm -rf $(LIBDIR) $(EXEDIR)

# END