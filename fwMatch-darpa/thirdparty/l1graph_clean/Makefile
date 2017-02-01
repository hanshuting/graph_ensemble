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
all: gibbs_mex.$(MEX_EXT)

# Handy variable reference (to implement Don't Repeat Yourself)
# $@ -- name of the target
# $< -- name of the first prerequisite
# $^ -- name of all prerequisites, separated by spaces

gibbs_mex.$(MEX_EXT): gibbs_mex.cpp
	mex -largeArrayDims -v -O $< -output $@

clean:
	rm -f *.$(MEX_EXT)

