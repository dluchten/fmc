# Standard config file for make, using gcc
#
# D. M. Luchtenburg
# Princeton University
#
CXX = g++

# Optimized settings
# CXXFLAGS = -O4 -funroll-loops -DNDEBUG

# Debugging
CXXFLAGS =  -g -Wall -Wextra -ansi -Woverloaded-virtual

lib_dirs = \
	-L/usr/local/lib

include_dirs = \
	-I/usr/local/include
