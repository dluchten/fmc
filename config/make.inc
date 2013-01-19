# Standard config file for make, using gcc
#
# D. M. Luchtenburg
# Princeton University
#
CXX = g++

# Optimized settings
# CXXFLAGS = -march=native -O2 -pipe -DNDEBUG
#CXXFLAGS = -march=core2 -msse4.2 -msse4.1 -mtune=core2 -O2 -pipe -DNDEBUG
CXXFLAGS = -O2 -pipe -DNDEBUG


# Debugging
#CXXFLAGS = -std=c++0x -g -Wall -Wextra -ansi -Woverloaded-virtual
#CXXFLAGS = -g -Wall -Wextra -ansi -Woverloaded-virtual

lib_dirs = \
	-L/usr/local/lib

include_dirs = \
	-I/usr/local/include
