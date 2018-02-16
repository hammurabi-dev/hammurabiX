# -*- Makefile -*-

default: hammurabi lib clean

hammurabi:
	$(MAKE) $(MFLAGS) -C install default

lib:
	$(MAKE) $(MFLAGS) -C install lib

doc:
	$(MAKE) $(MFLAGS) -C install doc

clean:
	$(MAKE) $(MFLAGS) -C install clean

wipe:
	$(MAKE) $(MFLAGS) -C install wipe

test: lib
	$(MAKE) $(MFLAGS) -C install clean
	$(MAKE) $(MFLAGS) -C unitest default
	$(MAKE) $(MFLAGS) -C unitest clean
	$(MAKE) $(MFLAGS) -C unitest run
	$(MAKE) $(MFLAGS) -C unitest wipe
	$(MAKE) $(MFLAGS) -C install wipe
