CXXFLAGS := -O3 -std=c++11
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

all: csaAssign.$(MEX_EXT) csaAssignPermMat.$(MEX_EXT) test

csaAssign.$(MEX_EXT): csa.cc csa.hh csaAssign.cc
	mex -largeArrayDims -v -g -O csaAssign.cc csa.cc -output csaAssign.$(MEX_EXT)

csaAssignPermMat.$(MEX_EXT): csa.cc csa.hh csaAssignPermMat.cc
	mex -I../../src/include -I/usr/local/include -L/usr/local/lib -largeArrayDims -v -g -O csaAssignPermMat.cc csa.cc -output csaAssignPermMat_mex.$(MEX_EXT)

test: csa.cc csa.hh test.cc
	$(CXX) csa.cc test.cc -o test

clean:
	rm test *.o *.$(MEX_EXT)
