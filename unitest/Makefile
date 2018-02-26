include ../install/makefile.ini

ifeq ($(CXX),)
CXX=g++
endif

ifeq ($(GSL),)
GSL=/usr
endif

ifeq ($(FFTW),)
FFTW=/usr/local
endif

ifeq ($(CFITSIO),)
CFITSIO=/usr/local
endif

ifeq ($(HEALPIX),)
$(error "missing HEALPIX")
endif

ifeq ($(HEALPIX_TARGET),)
$(error "missing HEALPIX_TARGET")
endif

HEALPIX_FULL = $(HEALPIX)/src/cxx/$(HEALPIX_TARGET)

SRCDIR=.
SOURCES=$(wildcard *.cc)
OBJDIR=./build
OBJ=$(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(SOURCES:.cc=.o))
LIBDIR=../lib
INCDIR=../include
HEADS=$(wildcard $(INCDIR)*.h)
EXEDIR=./bin
TESTEXE=$(addprefix $(EXEDIR)/,$(SOURCES:.cc=.x))

MRCDIR=../src
MBJDIR=../build
MOURCES=$(shell find $(MRCDIR) -type f -name "*.cc")
MBJ=$(patsubst $(MRCDIR)/%,$(MBJDIR)/%,$(MOURCES:.cc=.o))


# COMPILE OPTIONS
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -fopenmp -pthread -DDEBUG $(DEBUG)

EXT_INC = -I$(CFITSIO)/include -I$(GSL)/include -I$(FFTW)/include \
	-I$(HEALPIX_FULL)/include \
	-I $(INCDIR) -I.
CXXFLAGS += $(EXT_INC)
 
ALL_L = -L$(CFITSIO)/lib -L$(GSL)/lib -L$(FFTW)/lib \
	-L$(HEALPIX_FULL)/lib -L$(LIBDIR) -L.
LD_FLAGS = $(ALL_L) \
	-lhealpix_cxx -lcxxsupport -lc_utils -lfftpack -lcfitsio \
	-lgsl -lgslcblas -lm -lfftw3 -fopenmp $(DEBUG)
LFLAGS = $(ALL_L) -lhammurabi \
	-lhealpix_cxx -lcxxsupport -lc_utils -lfftpack -lcfitsio \
	-lgsl -lgslcblas -lm -lfftw3 -fopenmp $(DEBUG)

# RULES

default: $(LIBDIR)/libhammurabi.so $(TESTEXE)

$(LIBDIR)/libhammurabi.so: $(MBJ)
	mkdir -p $(LIBDIR)
	$(CXX) -o $@ -shared -fPIC $(LD_FLAGS) $^

$(MBJDIR)/%.o: $(MRCDIR)/%.cc $(HEADS)
	mkdir -p $(dir $@)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(TESTEXE): | $(EXEDIR)

$(EXEDIR):
	mkdir -p $(EXEDIR)

$(EXEDIR)/%.x: $(OBJDIR)/%.o
	$(CXX) $^ -o $@ $(LFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	mkdir -p $(dir $@)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

help:
	echo $(MBJ)

.PHONY: clean
clean:
	rm -rf $(OBJDIR) $(MBJDIR)

.PHONY: wipe
wipe:
	rm -rf $(EXEDIR) $(LIBDIR)

run: ./bin/param_unitest.x ./bin/int_unitest.x
	export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(LIBDIR)
	./bin/param_unitest.x
	./bin/int_unitest.x

#END
