all: run_tests

TEST_FILES= \
	main.o \
	legendrebasis_test.o \
	gpcexpansion_test.o \
	interp_trilin_test.o \
	velocityfield_snap_2d_test.o

include ../config/make.inc

CXXFLAGS += $(include_dirs) -I../src
LDFLAGS += $(lib_dirs)
BUILDDIR = ../build
MAINLIB = libmain.a
LIBS = $(BUILDDIR)/$(MAINLIB) -lgtest -lpthread -lm 

MAKEDEPEND = g++ -MM

.PHONY: clean distclean depend

run_tests: runner
	./runner 2> runner.err

runner: lib $(TEST_FILES)
	$(CXX) $(LDFLAGS) -o $@ $(TEST_FILES) $(LIBS)

lib:
	cd $(BUILDDIR) && make $(MAINLIB)

% : %.o
	$(CXX) $(LDFLAGS) -o $@ $< $(LIBS)

%.o : %.cc
	$(CXX) -c $(CXXFLAGS) $(INCPATH) $<

clean:
	@-$(RM) -rf *.o *~ *.err
	@-$(RM) runner

distclean: clean
	@-$(RM) runner

depend:
	$(MAKEDEPEND) $(CXXFLAGS) $(INCPATH) *.cc > .depend

-include .depend
