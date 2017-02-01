BOOST_DIR=../externalpackages/
UTILITIES_DIR=../utilities
DATA_DIR=../data
MODEL_DIR=../model
LOSS_DIR=../loss
SOLVER_DIR=../solver
INNER_SOLVER_DIR=../solver/bmrminnersolver

CORE_OBJ = sml.o common.o configuration.o bmrmexception.o timer.o
DATA_OBJ = data.o vecfeature.o veclabel.o vecdata.o varlenveclabel.o multilabelvecdata.o seqlabel.o seqfeature.o seqdata.o
MODEL_OBJ = model.o
LOSS_OBJ = loss.o wtamulticlassloss.o rocscoreloss.o scalarloss.o binaryclassificationloss.o ndcgrankloss.o lap.o rampndcgrankloss.o rankloss.o smmsegmentationloss.o #softmarginloss.o univariateregressionloss.o epsiloninsensitiveloss.o\
	    #preferencerankloss.o ordinalregressionloss.o rankloss.o

SOLVER_OBJ = rampbmrm.o
INNER_SOLVER_OBJ = bmrminnersolver.o l2n2_bmrmdualinnersolver.o l2n2_daifletcherpgm.o l2n2_prloqo.o l2n2_linesearch.o l2n2_gmn.o

FACTORY_SRC = ${LOSS_DIR}/ramplossfactory.hpp ${DATA_DIR}/datafactory.hpp ${INNER_SOLVER_DIR}/bmrminnersolverfactory.hpp

# CXX = g++

### for debugging ###
#CFLAGS = -Wall -O0 -g -DDEBUG -I${DATA_DIR} -I${MODEL_DIR} -I${UTILITIES_DIR} -I${LOSS_DIR} -I${SOLVER_DIR} -I${INNER_SOLVER_DIR} -I${BOOST_DIR}

### for release ###
CFLAGS = -Wall -O3 -DNDEBUG -I${DATA_DIR} -I${MODEL_DIR} -I${UTILITIES_DIR} -I${LOSS_DIR} -I${SOLVER_DIR} -I${INNER_SOLVER_DIR} -I${BOOST_DIR}

LDFLAGS = 


.PHONY: all

all: ramp-bmrm-train ramp-bmrm-predict

# executable
ramp-bmrm-train: ramp-bmrm-train.cpp ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ} ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} ramp-bmrm-train.cpp -o ramp-bmrm-train ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ} 

ramp-bmrm-predict: ramp-bmrm-predict.cpp ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} ramp-bmrm-predict.cpp -o ramp-bmrm-predict ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} 


# solver objects
rampbmrm.o: ${SOLVER_DIR}/rampbmrm.hpp ${SOLVER_DIR}/rampbmrm.cpp ${FACTORY_SRC}
	${CXX} ${CFLAGS} -c ${SOLVER_DIR}/rampbmrm.cpp

bmrminnersolver.o: ${INNER_SOLVER_DIR}/bmrminnersolver.hpp ${INNER_SOLVER_DIR}/bmrminnersolver.cpp
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/bmrminnersolver.cpp

l2n2_bmrmdualinnersolver.o: ${INNER_SOLVER_DIR}/l2n2_bmrmdualinnersolver.hpp ${INNER_SOLVER_DIR}/l2n2_bmrmdualinnersolver.cpp
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_bmrmdualinnersolver.cpp

l2n2_daifletcherpgm.o: ${INNER_SOLVER_DIR}/l2n2_daifletcherpgm.hpp ${INNER_SOLVER_DIR}/l2n2_daifletcherpgm.cpp
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_daifletcherpgm.cpp

l2n2_prloqo.o: ${INNER_SOLVER_DIR}/l2n2_prloqo.hpp ${INNER_SOLVER_DIR}/l2n2_prloqo.cpp
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_prloqo.cpp

l2n2_linesearch.o: ${INNER_SOLVER_DIR}/l2n2_linesearch.hpp ${INNER_SOLVER_DIR}/l2n2_linesearch.cpp
	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_linesearch.cpp

l2n2_gmn.o: ${INNER_SOLVER_DIR}/l2n2_gmn.hpp ${INNER_SOLVER_DIR}/l2n2_gmn.cpp
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

epsiloninsensitiveloss.o: ${LOSS_DIR}/epsiloninsensitiveloss.hpp ${LOSS_DIR}/epsiloninsensitiveloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/epsiloninsensitiveloss.cpp

wtamulticlassloss.o: ${LOSS_DIR}/wtamulticlassloss.hpp ${LOSS_DIR}/wtamulticlassloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/wtamulticlassloss.cpp

ndcgrankloss.o: ${LOSS_DIR}/ndcgrankloss.hpp ${LOSS_DIR}/ndcgrankloss.cpp lap.o
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/ndcgrankloss.cpp

rampndcgrankloss.o: ${LOSS_DIR}/rampndcgrankloss.hpp ${LOSS_DIR}/rampndcgrankloss.cpp
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/rampndcgrankloss.cpp

lap.o: ${LOSS_DIR}/lap.cpp ${LOSS_DIR}/lap.hpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/lap.cpp

preferencerankloss.o: ${LOSS_DIR}/preferencerankloss.hpp ${LOSS_DIR}/preferencerankloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/preferencerankloss.cpp

ordinalregressionloss.o: ${LOSS_DIR}/ordinalregressionloss.hpp ${LOSS_DIR}/ordinalregressionloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/ordinalregressionloss.cpp

rocscoreloss.o: ${LOSS_DIR}/rocscoreloss.hpp ${LOSS_DIR}/rocscoreloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/rocscoreloss.cpp

smmsegmentationloss.o: ${LOSS_DIR}/smmsegmentationloss.hpp ${LOSS_DIR}/smmsegmentationloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/smmsegmentationloss.cpp


# data objects
data.o: ${DATA_DIR}/data.hpp ${DATA_DIR}/data.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/data.cpp

vecfeature.o: ${DATA_DIR}/vecfeature.hpp ${DATA_DIR}/vecfeature.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/vecfeature.cpp

veclabel.o: ${DATA_DIR}/veclabel.hpp ${DATA_DIR}/veclabel.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/veclabel.cpp

varlenveclabel.o: ${DATA_DIR}/varlenveclabel.hpp ${DATA_DIR}/varlenveclabel.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/varlenveclabel.cpp

vecdata.o: ${DATA_DIR}/vecdata.hpp ${DATA_DIR}/vecdata.cpp veclabel.o vecfeature.o ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/vecdata.cpp

multilabelvecdata.o: ${DATA_DIR}/multilabelvecdata.hpp ${DATA_DIR}/multilabelvecdata.cpp varlenveclabel.o vecfeature.o ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/multilabelvecdata.cpp

seqlabel.o: ${DATA_DIR}/seqlabel.hpp ${DATA_DIR}/seqlabel.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/seqlabel.cpp

seqfeature.o: ${DATA_DIR}/seqfeature.hpp ${DATA_DIR}/seqfeature.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/seqfeature.cpp

seqdata.o: ${DATA_DIR}/seqdata.hpp ${DATA_DIR}/seqdata.cpp data.o seqlabel.o seqfeature.o ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/seqdata.cpp


# model objects
model.o: ${MODEL_DIR}/model.hpp ${MODEL_DIR}/model.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${MODEL_DIR}/model.cpp


# core objects
sml.o: ${UTILITIES_DIR}/sml.hpp ${UTILITIES_DIR}/sml.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/sml.cpp

#info.o: ${UTILITIES_DIR}/info.hpp ${UTILITIES_DIR}/info.cpp
#	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/info.cpp

common.o: ${UTILITIES_DIR}/common.hpp ${UTILITIES_DIR}/common.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/common.cpp

configuration.o: ${UTILITIES_DIR}/configuration.hpp ${UTILITIES_DIR}/configuration.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/configuration.cpp

bmrmexception.o: ${UTILITIES_DIR}/bmrmexception.hpp ${UTILITIES_DIR}/bmrmexception.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/bmrmexception.cpp

timer.o: ${UTILITIES_DIR}/timer.hpp ${UTILITIES_DIR}/timer.cpp
	${CXX} ${CFLAGS} -c ${UTILITIES_DIR}/timer.cpp

clean:
	rm *.o ramp-bmrm-train ramp-bmrm-predict;

clean-loss:
	rm *loss.o ramp-bmrm-train;

