BOOST_DIR=../externalpackages/
UTILITIES_DIR=../utilities
DATA_DIR=../data
MODEL_DIR=../model
LOSS_DIR=../loss
SOLVER_DIR=../solver
INNER_SOLVER_DIR= #../solver/bmrminnersolver

CORE_OBJ = sml.o common.o configuration.o bmrmexception.o timer.o
DATA_OBJ = data.o onlinevecdata.o vecdata.o
MODEL_OBJ = model.o
LOSS_OBJ = loss.o scalarloss.o binaryclassificationloss.o  softmarginloss.o 

SOLVER_OBJ = olbfgs.o lbfgs.o sublbfgs.o smd.o 
INNER_SOLVER_OBJ = #bmrminnersolver.o l2n2_bmrmdualinnersolver.o l2n2_daifletcherpgm.o l2n2_prloqo.o l2n2_linesearch.o # l2n2_smo.o l2n2_gmn.o

FACTORY_SRC = #${LOSS_DIR}/lossfactory.hpp ${DATA_DIR}/datafactory.hpp ${SOLVER_DIR}/solverfactory.hpp ${INNER_SOLVER_DIR}/bmrminnersolverfactory.hpp

# CXX = g++

### for debugging ###
#CFLAGS = -O0 -g -DDEBUG -I${DATA_DIR} -I${MODEL_DIR} -I${UTILITIES_DIR} -I${LOSS_DIR} -I${SOLVER_DIR} -I${INNER_SOLVER_DIR} -I${BOOST_DIR}

### for release ###
CFLAGS = -Wall -O3 -DNDEBUG -I${DATA_DIR} -I${MODEL_DIR} -I${UTILITIES_DIR} -I${LOSS_DIR} -I${SOLVER_DIR} -I${BOOST_DIR}

LDFLAGS = 


.PHONY: all

all: online-train sublbfgs-train
# executable
online-train: online-train.cpp ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ} ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} online-train.cpp -o online-train ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ} 


sublbfgs-train: sublbfgs-train.cpp ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ} ${FACTORY_SRC}
	${CXX} ${CFLAGS} ${LDFLAGS} sublbfgs-train.cpp -o sublbfgs-train ${CORE_OBJ} ${DATA_OBJ} ${MODEL_OBJ} ${LOSS_OBJ} ${SOLVER_OBJ} ${INNER_SOLVER_OBJ} 


# solver objects
sublbfgs.o: ${SOLVER_DIR}/sublbfgs.hpp ${SOLVER_DIR}/sublbfgs.cpp ${FACTORY_SRC}
	${CXX} ${CFLAGS} -c ${SOLVER_DIR}/sublbfgs.cpp

lbfgs.o: ${SOLVER_DIR}/lbfgs.hpp ${SOLVER_DIR}/lbfgs.cpp ${FACTORY_SRC}
	${CXX} ${CFLAGS} -c ${SOLVER_DIR}/lbfgs.cpp

olbfgs.o: ${SOLVER_DIR}/olbfgs.hpp ${SOLVER_DIR}/olbfgs.cpp ${FACTORY_SRC}
	${CXX} ${CFLAGS} -c ${SOLVER_DIR}/olbfgs.cpp

smd.o: ${SOLVER_DIR}/smd.hpp ${SOLVER_DIR}/smd.cpp ${FACTORY_SRC}
	${CXX} ${CFLAGS} -c ${SOLVER_DIR}/smd.cpp

# bmrm.o: ${SOLVER_DIR}/bmrm.hpp ${SOLVER_DIR}/bmrm.cpp ${FACTORY_SRC}
# 	${CXX} ${CFLAGS} -c ${SOLVER_DIR}/bmrm.cpp

# bmrminnersolver.o: ${INNER_SOLVER_DIR}/bmrminnersolver.hpp ${INNER_SOLVER_DIR}/bmrminnersolver.cpp
# 	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/bmrminnersolver.cpp

# l2n2_bmrmdualinnersolver.o: ${INNER_SOLVER_DIR}/l2n2_bmrmdualinnersolver.hpp ${INNER_SOLVER_DIR}/l2n2_bmrmdualinnersolver.cpp bmrminnersolver.o
# 	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_bmrmdualinnersolver.cpp

# l2n2_daifletcherpgm.o: ${INNER_SOLVER_DIR}/l2n2_daifletcherpgm.hpp ${INNER_SOLVER_DIR}/l2n2_daifletcherpgm.cpp l2n2_bmrmdualinnersolver.o
# 	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_daifletcherpgm.cpp

# l2n2_prloqo.o: ${INNER_SOLVER_DIR}/l2n2_prloqo.hpp ${INNER_SOLVER_DIR}/l2n2_prloqo.cpp l2n2_bmrmdualinnersolver.o
# 	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_prloqo.cpp

# l2n2_linesearch.o: ${INNER_SOLVER_DIR}/l2n2_linesearch.hpp ${INNER_SOLVER_DIR}/l2n2_linesearch.cpp
# 	${CXX} ${CFLAGS} -c ${INNER_SOLVER_DIR}/l2n2_linesearch.cpp



# loss objects
loss.o: ${LOSS_DIR}/loss.hpp ${LOSS_DIR}/loss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/loss.cpp

scalarloss.o: ${LOSS_DIR}/scalarloss.hpp ${LOSS_DIR}/scalarloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/scalarloss.cpp

binaryclassificationloss.o: ${LOSS_DIR}/binaryclassificationloss.hpp ${LOSS_DIR}/binaryclassificationloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/binaryclassificationloss.cpp

softmarginloss.o: ${LOSS_DIR}/softmarginloss.hpp ${LOSS_DIR}/softmarginloss.cpp 
	${CXX} ${CFLAGS} -c ${LOSS_DIR}/softmarginloss.cpp


# data objects
data.o: ${DATA_DIR}/data.hpp ${DATA_DIR}/data.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/data.cpp

onlinevecdata.o: ${DATA_DIR}/onlinevecdata.hpp ${DATA_DIR}/onlinevecdata.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/onlinevecdata.cpp

vecdata.o: ${DATA_DIR}/vecdata.hpp ${DATA_DIR}/vecdata.cpp ${CORE_OBJ}
	${CXX} ${CFLAGS} -c ${DATA_DIR}/vecdata.cpp


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
	rm *.o;

clean-loss:
	rm *loss.o;

