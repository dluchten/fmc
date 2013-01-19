# Makefile for main executable
#
# D. M. Luchtenburg
# Princeton University
#
MAINLIB = libmain.a
EXECUTABLES = fmc-gyre fmc-snaps

all: $(MAINLIB) $(EXECUTABLES)

include ../config/make.inc

OBJS = \
	stopwatch.o \
	integrator.o \
	heun.o \
	runge-kutta.o \
	doublegyre.o \
	orthopolybasis.o \
	legendrebasis.o \
	gpcexpansion.o \
	grid2d.o \
	flowmap2d.o \
	interp.o \
	interp_trilin.o \
	velocityfield_snap_2d.o

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

fmc-gyre : $(OBJS) fmc-gyre.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

fmc-snaps : $(OBJS) fmc-snaps.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

clean:
	@-$(RM) *.o
	@-$(RM) .depend
	@-$(RM) *~
	@-$(RM) $(EXECUTABLES) $(MAINLIB)
	@-$(RM) *.dat

distclean: clean
	@-$(RM) $(EXECUTABLES) $(MAINLIB)

depend:
	$(MAKEDEPEND) $(CXXFLAGS) ../src/*.cc > .depend

-include .depend
