# -*- Makefile -*-

default: hammurabi lib clean

hammurabi:
	$(MAKE) $(MFLAGS) -C INSTALL default

lib:
	$(MAKE) $(MFLAGS) -C INSTALL lib

doc:
	$(MAKE) $(MFLAGS) -C INSTALL doc

clean:
	$(MAKE) $(MFLAGS) -C INSTALL clean
	$(MAKE) $(MFLAGS) -C unitest clean

wipe:
	$(MAKE) $(MFLAGS) -C INSTALL wipe
	$(MAKE) $(MFLAGS) -C unitest wipe

test:
	$(MAKE) $(MFLAGS) -C unitest default
	$(MAKE) $(MFLAGS) -C unitest clean
	$(MAKE) $(MFLAGS) -C unitest run
	$(MAKE) $(MFLAGS) -C unitest wipe

install:
	$(MAKE) $(MFLAGS) -C INSTALL install
	$(MAKE) $(MFLAGS) -C INSTALL wipe
	$(MAKE) $(MFLAGS) -C unitest wipe