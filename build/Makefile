# Makefile for main executable
#
# D. M. Luchtenburg
# Princeton University
#
MAINLIB = libmain.a
EXECUTABLES = main

all: $(MAINLIB) $(EXECUTABLES)

include ../config/make.inc

OBJS = \
	main.o \
	heun.o \
	doublegyre.o \
	gpcbasis.o \
	legendrebasis.o

VPATH = ../src

AR = ar
ARFLAGS = -r

MAKEDEPEND = g++ -MM

LDLIBS = -lm
LDFLAGS += $(lib_dirs)
CXXFLAGS += $(include_dirs)

.PHONY: clean distclean depend

$(MAINLIB): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

$(EXECUTABLES) : $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	@-$(RM) *.o
	@-$(RM) .depend
	@-$(RM) *~
	@-$(RM) $(EXECUTABLES) $(MAINLIB)

distclean: clean
	@-$(RM) $(EXECUTABLES) $(MAINLIB)

depend:
	$(MAKEDEPEND) $(CXXFLAGS) ../src/*.cc > .depend

-include .depend
