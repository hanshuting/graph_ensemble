DIRS := . MinCost GEOM 

SOURCES := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cpp))
OBJS := $(patsubst %.cpp, %.o, $(SOURCES))

CFLAGS := -DPERFECT_MATCHING_DOUBLE -O3 -D_NDEBUG -fPIC
#CFLAGS := -DPERFECT_MATCHING_DOUBLE -O0 -g

UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
	MEX_EXT := mexa64
	CC      := g++
endif
ifeq ($(UNAME), Darwin)
	MEX_EXT := mexmaci64
	CC      := clang++
endif

#LIBS :=  -lc++
INCLUDES := 
LIBDIR := 

all: blossom_double_mex.$(MEX_EXT) blossom_double_mex_debug.$(MEX_EXT) blossom_sparse_mex.$(MEX_EXT) blossom_sparse_mex_debug.$(MEX_EXT)

blossom_sparse_mex.$(MEX_EXT): blossom_sparse_mex.cc mexcpp.h ${OBJS}
	LD=$(CC) mex -DPERFECT_MATCHING_DOUBLE -largeArrayDims -v -O blossom_sparse_mex.cc -output blossom_sparse_mex.$(MEX_EXT) ${OBJS}

blossom_sparse_mex_debug.$(MEX_EXT): blossom_sparse_mex.cc mexcpp.h ${OBJS}
	LD=$(CC) mex -DPERFECT_MATCHING_DOUBLE -largeArrayDims -v -g blossom_sparse_mex.cc -output blossom_sparse_mex_debug.$(MEX_EXT) ${OBJS}

blossom_double_mex.$(MEX_EXT): blossom_double_mex.cc mexcpp.h ${OBJS}
	LD=$(CC) mex -DPERFECT_MATCHING_DOUBLE -largeArrayDims -v -O blossom_double_mex.cc -output blossom_double_mex.$(MEX_EXT) ${OBJS}

blossom_double_mex_debug.$(MEX_EXT): blossom_double_mex.cc mexcpp.h ${OBJS}
	LD=$(CC) mex -DPERFECT_MATCHING_DOUBLE -largeArrayDims -v -g blossom_double_mex.cc -output blossom_double_mex_debug.$(MEX_EXT) ${OBJS}

.cpp.o:
	$(CC) $(CFLAGS) ${INCLUDES} $< -c -o $@

clean:
	rm -f ${OBJS} blossom5 *.$(MEX_EXT)
