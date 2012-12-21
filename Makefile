# Main makefile
#
# D. M. Luchtenburg
# Princeton University
#
DIRS = build test doc

main:
	cd build && $(MAKE)

test:
	cd test && $(MAKE)

doc:
	cd doc && $(MAKE)

all: main test doc

clean:
	@$(RM) *~
	@for dir in $(DIRS); do ( cd $$dir && $(MAKE) clean; ) done

distclean: clean
	@for dir in $(DIRS); do \
	  ( cd $$dir && $(MAKE) distclean; )\
	done

.PHONY: main test doc clean distclean

