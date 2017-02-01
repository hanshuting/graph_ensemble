CXX      := clang++
CXXFLAGS := -O3 -g -std=c++11 -ffast-math -DNDEBUG
INCFLAGS := -I../../thirdparty -I../include -I/usr/local/include
#CXXFLAGS := -O3 -g -std=c++11 -ffast-math -framework Accelerate
#CXXFLAGS := -O3 -g -DNDEBUG -std=c++11

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	MEX_EXT := mexa64
endif
ifeq ($(UNAME), Darwin)
	MEX_EXT := mexmaci64
endif

all: fwBipartite_mex.$(MEX_EXT) fwBipartite_mex_debug.$(MEX_EXT)

.PHONY: clean

clean:
	rm *.o *.$(MEX_EXT)

# TODO: Wrap the dlibLap only.

fwBipartite_mex.$(MEX_EXT): fwBipartite_mex.cpp
	mex $(INCFLAGS) -v -largeArrayDims -O -g fwBipartite_mex.cpp -o fwBipartite_mex.$(MEX_EXT)

fwBipartite_mex_debug.$(MEX_EXT): fwBipartite_mex.cpp 
	mex $(INCFLAGS) -v -largeArrayDims -g fwBipartite_mex.cpp -o fwBipartite_mex_debug.$(MEX_EXT)

