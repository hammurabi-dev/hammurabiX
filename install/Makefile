include ./makefile.ini

SRCDIR=../src
OBJDIR=../build
SOURCES=$(shell find $(SRCDIR) -type f -name *.cc)
OBJ=$(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(SOURCES:.cc=.o))
INCDIR=../include
HEADS=$(wildcard $(INCDIR)*.h)
LIBDIR=../lib
EXEDIR=../bin

# COMPILE OPTIONS
CXXFLAGS = -O3 -Wall -Wextra -pedantic -fPIC -Wno-deprecated -fopenmp -pthread $(DEBUG)

EXT_INC = -I$(CFITSIO)/include -I$(GSL)/include -I$(FFTW)/include \
	-I$(HEALPIX)/include \
	-I $(INCDIR) -I.
CXXFLAGS += $(EXT_INC)
 
ALL_L = -L$(CFITSIO)/lib -L$(GSL)/lib -L$(FFTW)/lib \
	-L$(HEALPIX)/lib\
	-L$(LIBDIR) -L.
LD_FLAGS = $(ALL_L) \
	-lhealpix_cxx -lcxxsupport -lc_utils -lfftpack -lcfitsio \
	-lgsl -lgslcblas -lm -lfftw3 -fopenmp
LFLAGS = $(LD_FLAGS)

# RULES

default: $(EXE)

$(EXE): | $(EXEDIR)

$(EXEDIR):
	mkdir -p $(EXEDIR)

$(EXE): $(OBJ)
	$(CXX) $^ -o $(EXEDIR)/$(EXE) $(LFLAGS)
	cp ../templates/* $(EXEDIR)

$(OBJDIR)%.o: $(SRCDIR)%.cc $(HEADS)
	mkdir -p $(dir $@)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(LIBDIR)/libhammurabi.so: $(OBJ)
	mkdir -p $(LIBDIR)
	$(CXX) -o $@ -shared $(LD_FLAGS) $^

help:
	echo $(OBJ)

documentation: install/doxy.in
	doxygen $^

dynamiclib: $(LIBDIR)/libhammurabi.so

.PHONY: clean
clean:
	rm -rf $(LIBDIR); rm -rf $(OBJDIR); rm -rf $(EXEDIR); rm -rf html latex;

.PHONY: tar 
tar:
	tar cfv HamX.tar.gz AUTHROS COPYING readme.md ../templates $(INCDIR) $(OBJDIR) $(LIBDIR)
# END
