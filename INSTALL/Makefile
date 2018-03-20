include ./makefile.ini

ifeq ($(INS_PATH),)
INS_PATH=/usr/local/hammurabi
endif

ifeq ($(CXX),)
CXX=g++
endif

ifeq ($(DEBUG),)
DEBUG= -DNDEBUG
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

SRCDIR=../src
OBJDIR=../build
SOURCES=$(shell find $(SRCDIR) -type f -name "*.cc")
OBJ=$(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(SOURCES:.cc=.o))
INCDIR=../include
HEADS=$(wildcard $(INCDIR)*.h)
LIBDIR=../lib
EXEDIR=../bin

# COMPILE OPTIONS
CXXFLAGS = -std=c++11 -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -fopenmp -pthread $(DEBUG)

EXT_INC = -I$(CFITSIO)/include -I$(GSL)/include -I$(FFTW)/include \
	-I$(HEALPIX_FULL)/include \
	-I $(INCDIR) -I.
CXXFLAGS += $(EXT_INC)
 
ALL_L = -L$(CFITSIO)/lib -L$(GSL)/lib -L$(FFTW)/lib \
	-L$(HEALPIX_FULL)/lib
LD_FLAGS = $(ALL_L) \
	-lhealpix_cxx -lcxxsupport -lc_utils -lfftpack -lcfitsio \
	-lgsl -lgslcblas -lfftw3_omp -lfftw3 -lm -fopenmp $(DEBUG)
LFLAGS = $(LD_FLAGS)

# RULES

default: $(EXE)

$(EXE): | $(EXEDIR)

$(EXEDIR):
	mkdir -p $(EXEDIR)

$(EXE): $(OBJ)
	$(CXX) $^ -o $(EXEDIR)/$(EXE) $(LFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc $(HEADS)
	mkdir -p $(dir $@)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(LIBDIR)/libhammurabi.so: $(OBJ)
	mkdir -p $(LIBDIR)
	$(CXX) -o $@ -shared -fPIC $(LD_FLAGS) $^

.PHONY: help
help:
	echo $(OBJ)

.PHONY: doc
doc: install/doxy.in
	doxygen $^

.PHONY: lib
lib: $(LIBDIR)/libhammurabi.so


.PHONY: clean
clean:
	rm -rf $(OBJDIR)

.PHONY: wipe
wipe:
	rm -rf $(LIBDIR) $(EXEDIR) thml latex

.PHONY: install
install:
	mkdir -p $(INS_PATH)/bin
	mkdir -p $(INS_PATH)/lib
	mkdir -p $(INS_PATH)/include
	cp $(EXEDIR)/* $(INS_PATH)/bin; cp $(LIBDIR)/* $(INS_PATH)/lib; cp $(INCDIR)/* $(INS_PATH)/include

# END
