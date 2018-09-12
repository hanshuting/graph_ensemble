
UTILITIES_DIR=../../utilities
SML_DIR=../../linear-bmrm
DATACONTAINER_DIR=../../datacontainer

CORE_OBJ = sml.o info.o common.o configuration.o bmrmexception.o
DATA_OBJ = datacontainer.o veclabelvecfeaturedc.o

COMPILER = g++
CFLAGS = -O3 -DNDEBUG -Wall -I${DATACONTAINER_DIR} -I${UTILITIES_DIR} -I${SML_DIR}




test_veclabelvecfeaturedc: test_veclabelvecfeaturedc.cpp  ${CORE_OBJ} ${DATA_OBJ}
	${COMPILER} ${CFLAGS} test_veclabelvecfeaturedc.cpp -o test_veclabelvecfeaturedc ${CORE_OBJ} ${DATA_OBJ}


# data objects
veclabelvecfeaturedc.o: ${DATACONTAINER_DIR}/veclabelvecfeaturedc.hpp ${DATACONTAINER_DIR}/veclabelvecfeaturedc.cpp ${CORE_OBJ}
	${COMPILER} ${CFLAGS} -c ${DATACONTAINER_DIR}/veclabelvecfeaturedc.cpp

datacontainer.o: ${DATACONTAINER_DIR}/datacontainer.hpp ${DATACONTAINER_DIR}/datacontainer.cpp ${CORE_OBJ}
	${COMPILER} ${CFLAGS} -c ${DATACONTAINER_DIR}/datacontainer.cpp


# core objects
sml.o: ${SML_DIR}/sml.hpp ${SML_DIR}/sml.cpp
	${COMPILER} ${CFLAGS} -c ${SML_DIR}/sml.cpp

info.o: ${UTILITIES_DIR}/info.hpp ${UTILITIES_DIR}/info.cpp
	${COMPILER} ${CFLAGS} -c ${UTILITIES_DIR}/info.cpp

common.o: ${UTILITIES_DIR}/common.hpp ${UTILITIES_DIR}/common.cpp
	${COMPILER} ${CFLAGS} -c ${UTILITIES_DIR}/common.cpp

configuration.o: ${UTILITIES_DIR}/configuration.hpp ${UTILITIES_DIR}/configuration.cpp
	${COMPILER} ${CFLAGS} -c ${UTILITIES_DIR}/configuration.cpp

bmrmexception.o: ${UTILITIES_DIR}/bmrmexception.hpp ${UTILITIES_DIR}/bmrmexception.cpp
	${COMPILER} ${CFLAGS} -c ${UTILITIES_DIR}/bmrmexception.cpp


clean:
	rm *.o test_veclabelvecfeaturedc; mv *~ ~/tmp/


