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

SOURCES=$(wildcard *.cc)
OBJDIR=.
OBJ=$(SOURCES:.cc=.o)
LIBDIR=../lib
INCDIR=../include
EXEDIR=./bin
TESTEXE=$(addprefix $(EXEDIR)/,$(SOURCES:.cc=.x))

# COMPILE OPTIONS
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -fopenmp -pthread -g -pg -DDEBUG

EXT_INC = -I$(CFITSIO)/include -I$(GSL)/include -I$(FFTW)/include \
	-I$(HEALPIX_FULL)/include \
	-I $(INCDIR) -I.
CXXFLAGS += $(EXT_INC)
 
ALL_L = -L$(CFITSIO)/lib -L$(GSL)/lib -L$(FFTW)/lib \
	-L$(HEALPIX_FULL)/lib -L$(LIBDIR) -L.
LD_FLAGS = $(ALL_L) \
	-lhealpix_cxx -lcxxsupport -lc_utils -lfftpack -lcfitsio \
	-lgsl -lgslcblas -lm -lfftw3 -fopenmp -lhammurabi $(DEBUG)
LFLAGS = $(LD_FLAGS)

# RULES

default: $(TESTEXE)

$(TESTEXE): | $(EXEDIR)

$(EXEDIR):
	mkdir -p $(EXEDIR)

$(EXEDIR)/%.x: $(OBJDIR)/%.o
	$(CXX) $^ -o $@ $(LFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	$(CXX) -c $< -o $@ $(CXXFLAGS)

help:
	echo $(OBJ)

.PHONY: clean
clean:
	rm -rf $(OBJDIR)/*.o

.PHONY: wipe
wipe:
	rm -rf $(EXEDIR)

run: ./bin/param_unitest.x
	./bin/param_unitest.x

#END
