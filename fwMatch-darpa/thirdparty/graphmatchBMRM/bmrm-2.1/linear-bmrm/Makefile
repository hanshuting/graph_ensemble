BOOST_DIR=../externalpackages/
UTILITIES_DIR=../utilities
DATA_DIR=../data
MODEL_DIR=../model
LOSS_DIR=../loss
SOLVER_DIR=../solver
INNER_SOLVER_DIR=../solver/bmrminnersolver

ASP_DIR=asp

CORE_OBJ = sml.o common.o configuration.o bmrmexception.o timer.o
DATA_OBJ = data.o vecfeature.o veclabel.o varlenveclabel.o vecdata.o onlinevecdata.o multilabelvecdata.o seqlabel.o seqfeature.o seqdata.o graphdata.o
MODEL_OBJ = model.o
LOSS_OBJ = loss.o scalarloss.o binaryclassificationloss.o softmarginloss.o squaredsoftmarginloss.o rocscoreloss.o univariateregressionloss.o epsiloninsensitiveloss.o\
		   wtamulticlassloss.o expmulticlassloss.o preferencerankloss.o ordinalregressionloss.o rankloss.o multilabelloss.o\
	           smmsegmentationloss.o faceloss.o graphmatchloss.o

SOLVER_OBJ = bmrm.o
INNER_SOLVER_OBJ = bmrminnersolver.o l2n2_bmrmdualinnersolver.o l2n2_daifletcherpgm.o l2n2_prloqo.o l2n2_linesearch.o l2n2_gmn.o # l2n2_smo.o 

FACTORY_SRC = ${LOSS_DIR}/lossfactory.hpp ${DATA_DIR}/datafactory.hpp ${SOLVER_DIR}/solverfactory.hpp ${INNER_SOLVER_DIR}/bmrminnersolverfactory.hpp

X11_INC_DIR = /opt/X11/include/
X11_LIB_DIR = /opt/X11/lib/

# Use actual gcc (http://hpc.sourceforge.net/) instead of clang.
# CXX = g++
# CC = gcc

### for debugging ###

#CFLAGS = -Wall -O0 -g -DDEBUG -I${DATA_DIR} -I${MODEL_DIR} -I${UTILITIES_DIR} -I${LOSS_DIR} -I${SOLVER_DIR} -I${INNER_SOLVER_DIR} -I${BOOST_DIR} -I${X11_INC_DIR} -I${ASP_DIR} -L${X11_LIB_DIR}

### for release ###
CFLAGS = -Wall -O3 -DNDEBUG -I. -I${DATA_DIR} -I${MODEL_DIR} -I${UTILITIES_DIR} -I${LOSS_DIR} -I${SOLVER_DIR} -I${INNER_SOLVER_DIR} -I${BOOST_DIR} -I${X11_INC_DIR} -L${X11_LIB_DIR} -I${ASP_DIR} #-lgomp
#CFLAGS = -Wall -O3 -DNDEBUG -I${DATA_DIR} -I${MODEL_DIR} -I${UTILITIES_DIR} -I${LOSS_DIR} -I${SOLVER_DIR} -I${INNER_SOLVER_DIR} -I${BOOST_DIR} -I${ASP_DIR} -lgomp

LDFLAGS = 


.PHONY: all

all: linear-bmrm-train linear-bmrm-predict

# executable
linear-bmrm-train: linear-bmrm-train.cpp ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ} ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} linear-bmrm-train.cpp -o linear-bmrm-train ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ} lap.o gap.o -lm -lpthread -lX11

online-train: online-train.cpp ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ} ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} online-train.cpp -o online-train ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ} 

linear-bmrm-predict: linear-bmrm-predict.cpp ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} linear-bmrm-predict.cpp -o linear-bmrm-predict ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} lap.o gap.o -lm -lpthread -lX11


# solver objects
bmrm.o: ${SOLVER_DIR}/bmrm.hpp ${SOLVER_DIR}/bmrm.cpp ${FACTORY_SRC}
	${CXX} ${CFLAGS} -c ${SOLVER_DIR}/bmrm.cpp

bmrminnersolver.o: ${INNER_SOLVER_DIR}/bmrminnersolver.hpp ${INNER_SOLVER_DIR}/bmrminnersolver.cpp
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/bmrminnersolver.cpp

l2n2_bmrmdualinnersolver.o: ${INNER_SOLVER_DIR}/l2n2_bmrmdualinnersolver.hpp ${INNER_SOLVER_DIR}/l2n2_bmrmdualinnersolver.cpp bmrminnersolver.o
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_bmrmdualinnersolver.cpp

l2n2_daifletcherpgm.o: ${INNER_SOLVER_DIR}/l2n2_daifletcherpgm.hpp ${INNER_SOLVER_DIR}/l2n2_daifletcherpgm.cpp l2n2_bmrmdualinnersolver.o
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_daifletcherpgm.cpp

l2n2_prloqo.o: ${INNER_SOLVER_DIR}/l2n2_prloqo.hpp ${INNER_SOLVER_DIR}/l2n2_prloqo.cpp l2n2_bmrmdualinnersolver.o
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_prloqo.cpp

l2n2_linesearch.o: ${INNER_SOLVER_DIR}/l2n2_linesearch.hpp ${INNER_SOLVER_DIR}/l2n2_linesearch.cpp
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_linesearch.cpp

l2n2_gmn.o: ${INNER_SOLVER_DIR}/l2n2_gmn.hpp ${INNER_SOLVER_DIR}/l2n2_gmn.cpp l2n2_bmrmdualinnersolver.o
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_gmn.cpp


# loss objects
loss.o: ${LOSS_DIR}/loss.hpp ${LOSS_DIR}/loss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/loss.cpp

scalarloss.o: ${LOSS_DIR}/scalarloss.hpp ${LOSS_DIR}/scalarloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/scalarloss.cpp

binaryclassificationloss.o: ${LOSS_DIR}/binaryclassificationloss.hpp ${LOSS_DIR}/binaryclassificationloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/binaryclassificationloss.cpp

rankloss.o: ${LOSS_DIR}/rankloss.hpp ${LOSS_DIR}/rankloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/rankloss.cpp

univariateregressionloss.o: ${LOSS_DIR}/univariateregressionloss.hpp ${LOSS_DIR}/univariateregressionloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/univariateregressionloss.cpp

softmarginloss.o: ${LOSS_DIR}/softmarginloss.hpp ${LOSS_DIR}/softmarginloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/softmarginloss.cpp

squaredsoftmarginloss.o: ${LOSS_DIR}/squaredsoftmarginloss.hpp ${LOSS_DIR}/squaredsoftmarginloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/squaredsoftmarginloss.cpp

rocscoreloss.o: ${LOSS_DIR}/rocscoreloss.hpp ${LOSS_DIR}/rocscoreloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/rocscoreloss.cpp

epsiloninsensitiveloss.o: ${LOSS_DIR}/epsiloninsensitiveloss.hpp ${LOSS_DIR}/epsiloninsensitiveloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/epsiloninsensitiveloss.cpp

wtamulticlassloss.o: ${LOSS_DIR}/wtamulticlassloss.hpp ${LOSS_DIR}/wtamulticlassloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/wtamulticlassloss.cpp

expmulticlassloss.o: ${LOSS_DIR}/expmulticlassloss.hpp ${LOSS_DIR}/expmulticlassloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/expmulticlassloss.cpp

lap.o: ${LOSS_DIR}/lap.cpp ${LOSS_DIR}/lap.hpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/lap.cpp

preferencerankloss.o: ${LOSS_DIR}/preferencerankloss.hpp ${LOSS_DIR}/preferencerankloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/preferencerankloss.cpp

ordinalregressionloss.o: ${LOSS_DIR}/ordinalregressionloss.hpp ${LOSS_DIR}/ordinalregressionloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/ordinalregressionloss.cpp

faceloss.o: ${LOSS_DIR}/faceloss.cpp ${LOSS_DIR}/faceloss.hpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/faceloss.cpp

multilabelloss.o: ${LOSS_DIR}/multilabelloss.cpp ${LOSS_DIR}/multilabelloss.hpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/multilabelloss.cpp

smmsegmentationloss.o: ${LOSS_DIR}/smmsegmentationloss.cpp ${LOSS_DIR}/smmsegmentationloss.hpp seqdata.o
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/smmsegmentationloss.cpp

graphmatchloss.o: ${LOSS_DIR}/graphmatchloss.hpp ${LOSS_DIR}/graphmatchloss.cpp ${DATA_DIR}/graphdata.cpp gap.o lap.o ${DATA_DIR}/graphdata.cpp ${LOSS_DIR}/lap.cpp ${ASP_DIR}/gap.c
	${CXX} -fopenmp ${CFLAGS} -c ${LOSS_DIR}/graphmatchloss.cpp ${LOSS_DIR}/lap.cpp ${ASP_DIR}/gap.c



# data objects
data.o: ${DATA_DIR}/data.hpp ${DATA_DIR}/data.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/data.cpp

vecfeature.o: ${DATA_DIR}/vecfeature.hpp ${DATA_DIR}/vecfeature.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/vecfeature.cpp

veclabel.o: ${DATA_DIR}/veclabel.hpp ${DATA_DIR}/veclabel.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/veclabel.cpp

varlenveclabel.o: ${DATA_DIR}/varlenveclabel.hpp ${DATA_DIR}/varlenveclabel.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/varlenveclabel.cpp
	
onlinevecdata.o: ${DATA_DIR}/onlinevecdata.hpp ${DATA_DIR}/onlinevecdata.cpp vecdata.o ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/onlinevecdata.cpp

vecdata.o: ${DATA_DIR}/vecdata.hpp ${DATA_DIR}/vecdata.cpp data.o vecfeature.o veclabel.o ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/vecdata.cpp

multilabelvecdata.o: ${DATA_DIR}/multilabelvecdata.hpp ${DATA_DIR}/multilabelvecdata.cpp data.o vecfeature.o varlenveclabel.o ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/multilabelvecdata.cpp
	
	
seqlabel.o: ${DATA_DIR}/seqlabel.hpp ${DATA_DIR}/seqlabel.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/seqlabel.cpp

seqfeature.o: ${DATA_DIR}/seqfeature.hpp ${DATA_DIR}/seqfeature.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/seqfeature.cpp

seqdata.o: ${DATA_DIR}/seqdata.hpp ${DATA_DIR}/seqdata.cpp data.o seqlabel.o seqfeature.o ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/seqdata.cpp

graphdata.o: ${DATA_DIR}/graphdata.hpp ${DATA_DIR}/graphdata.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/graphdata.cpp



# model objects
model.o: ${MODEL_DIR}/model.hpp ${MODEL_DIR}/model.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${MODEL_DIR}/model.cpp


# Assignment solvers.

gap.o: ${ASP_DIR}/gap.h ${ASP_DIR}/gap.c
	${CXX} ${CFLAGS} -c ${ASP_DIR}/gap.c



# core objects
sml.o: ${UTILITIES_DIR}/sml.hpp ${UTILITIES_DIR}/sml.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/sml.cpp

common.o: ${UTILITIES_DIR}/common.hpp ${UTILITIES_DIR}/common.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/common.cpp

configuration.o: ${UTILITIES_DIR}/configuration.hpp ${UTILITIES_DIR}/configuration.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/configuration.cpp

bmrmexception.o: ${UTILITIES_DIR}/bmrmexception.hpp ${UTILITIES_DIR}/bmrmexception.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/bmrmexception.cpp

timer.o: ${UTILITIES_DIR}/timer.hpp ${UTILITIES_DIR}/timer.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/timer.cpp

clean:
	rm *.o linear-bmrm-train linear-bmrm-predict;

clean-loss:
	rm *loss.o linear-bmrm-train;

