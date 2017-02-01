CXX := /usr/local/bin/g++
CXXFLAGS := -v -O3
#CXXFLAGS := -O3 -g -std=c++11 -ffast-math -framework Accelerate
#CXXFLAGS := -O3 -g -DNDEBUG -std=c++11
SRCS := MRFEnergy.cpp minimize.cpp ordering.cpp treeProbabilities.cpp example.cpp
OBJS := $(patsubst %.cpp, %.o, $(SRCS))

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	MEX_EXT := mexa64
	CXX     := g++
endif
ifeq ($(UNAME), Darwin)
	MEX_EXT := mexmaci64
	CXX     := clang++
endif

all: example

example: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

#csaAssign.$(MEX_EXT): csa.cc csa.hh csaAssign.cc
#	mex -largeArrayDims -v -g -O csaAssign.cc csa.cc -o csaAssign.$(MEX_EXT)
#
#csaAssignPermMat.$(MEX_EXT): csa.cc csa.hh csaAssignPermMat.cc
#	mex -I../../src/include -I/usr/local/include -L/usr/local/lib -largeArrayDims -v -g -O csaAssignPermMat.cc csa.cc -o csaAssignPermMat_mex.$(MEX_EXT)
#
clean:
	rm example *.o *.$(MEX_EXT)
