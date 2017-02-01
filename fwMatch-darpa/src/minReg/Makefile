CXXFLAGS := -O3 -g
INCFLAGS := -DHAVE_EIGEN -I/usr/local/include -I/usr/local/include/eigen3 -I../../thirdparty/csa++ -I../include
#LDFLAGS  := -lboost_log_setup-mt -lboost_log-mt
LDFLAGS  := -L/usr/local/lib -lboost_log_setup-mt -lboost_log-mt
SOURCES  := fwMinRegFeats_mex_private.cc ../../thirdparty/csa++/csa.cc
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

all: fwMinRegFeats_mex.$(MEX_EXT) fwMinRegFeats_mex_debug.$(MEX_EXT)

fwMinRegFeats_mex.$(MEX_EXT): $(SOURCES)
	mex $(INCFLAGS) $(LDFLAGS) -largeArrayDims -v -O -g $(SOURCES) -o fwMinRegFeats_mex_private.$(MEX_EXT)

fwMinRegFeats_mex_debug.$(MEX_EXT): $(SOURCES)
	mex $(INCFLAGS) $(LDFLAGS) -largeArrayDims -v -g    $(SOURCES) -o fwMinRegFeats_mex_private_debug.$(MEX_EXT)

.PHONY: clean

clean:
	rm test *.o *.$(MEX_EXT)
