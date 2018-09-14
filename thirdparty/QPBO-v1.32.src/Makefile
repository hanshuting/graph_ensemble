#CXXFLAGS := -O3
CXXFLAGS := -ansi -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -pthread -std=c++0x
INCFLAGS := -I../../src/include -I/usr/local/include -L/usr/local/lib 

#CXXFLAGS := -O3 -g -std=c++11 -ffast-math -framework Accelerate
#CXXFLAGS := -O3 -g -DNDEBUG -std=c++11

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	MEX_EXT := mexa64
	CXX     := g++
endif
ifeq ($(UNAME), Darwin)
	MEX_EXT := mexmaci64
	CXX     := clang++
endif

# replace the following variable with your matlab mex path
MEX = /Applications/MATLAB_R2016a.app/bin/mex

all: example QPBO_double_mex.$(MEX_EXT) QPBO_double_mex_debug.$(MEX_EXT) QPBO_single_mex.$(MEX_EXT) QPBO_single_mex_debug.$(MEX_EXT)

example: example.cpp
	$(CXX) -v -o $@ example.cpp 

QPBO_double_mex.$(MEX_EXT): QPBO_mex.cpp
	$(MEX) $(INCFLAGS) CXXFLAGS='$(CXXFLAGS)' -DCOST_TYPE=double -largeArrayDims -v -O -g $? -output $@
#	mex $(INCFLAGS) CXXFLAGS='$(CXXFLAGS)' -DCOST_TYPE=double -largeArrayDims -v -O -g $? -output $@

QPBO_double_mex_debug.$(MEX_EXT): QPBO_mex.cpp
	$(MEX) $(INCFLAGS) CXXFLAGS='$(CXXFLAGS)' -DCOST_TYPE=double -largeArrayDims -v -g $? -output $@

QPBO_single_mex.$(MEX_EXT): QPBO_mex.cpp
	$(MEX) $(INCFLAGS) CXXFLAGS='$(CXXFLAGS)' -DCOST_TYPE=float -largeArrayDims -v -O -g $? -output $@

QPBO_single_mex_debug.$(MEX_EXT): QPBO_mex.cpp
	$(MEX) $(INCFLAGS) CXXFLAGS='$(CXXFLAGS)' -DCOST_TYPE=float -largeArrayDims -v -g $? -output $@

clean:
	rm -rf example *.$(MEX_EXT) *.dSYM/ 
