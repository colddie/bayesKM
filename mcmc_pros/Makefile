
# core compiling options
CXX = g++

CXX_STD = -std=c++11
OPT_FLAGS = -Wall -O3 -march=native -ffp-contract=fast -DARMA_NO_DEBUG

ARMA_INCLUDE_PATH = /home/tsun/bin/armadillo-9.400.4/install/include

# install location
INSTALL_PATH=/home/tsun/bin/mcmc-master/install/lib

# source directories
MCMC_HEADER_DIR = /home/tsun/bin/mcmc-master/install/include
MCMC_TEST_DIR = .

# general flags
LDFLAGS = -shared
CXXFLAGS = -fPIC -fopenmp $(CXX_STD) $(OPT_FLAGS) -I$(ARMA_INCLUDE_PATH) -I$(MCMC_HEADER_DIR)
LIBS= -L/home/tsun/bin/mcmc-master/install/lib -ldl -lmcmc /cm/shared/apps/openblas/0.2.20/lib/libopenblas.so


# Optim Test Files
# SOURCES_MCMC := rwmh_tac.cpp  # $(shell find $(MCMC_TEST_DIR) -name '*.cpp')

SRC = $(wildcard *.cpp)
OBJ = $(src:.cxx=.o)

TARGET = rwmh_tac

# $(TARGET): $(OBJ)
# 	$(CXX) $(CXXFLAGS) -o  $@ $^ $(SRC) $(LIBS) 
	
#sim1cm.o: sim1cm.c tpccm.h 
#sim2cm.o: sim2cm.c tpccm.h 
#simrtcm.o: simrtcm.c tpccm.h 
#simpct.o:
# $(TARGET).o: rwmh_tac_1tpc.cpp rwmh_tac_2tpc.cpp
#OBJECTS_MCMC := rwmh_tac_1tpc.o rwmh_tac_2tpc.o sim1cm.o sim2cm.o simrtcm.o #$(addprefix $(MCMC_TEST_DIR)/,$(SOURCES_MCMC:%.cpp=%.test))

#all: lib$(TARGET).so lib$(TARGET).a  #$(TARGET)  #$(addprefix $(MCMC_TEST_DIR)/,$(SOURCES_MCMC:%.cpp=%.test))    ##

# # core Optim files
# $(MCMC_TEST_DIR)/%.test: $(MCMC_TEST_DIR)/%.cpp 
# 	$(CXX) $(CXXFLAGS) $< -o $@ $(LIBS)

# models in a library                 ##
# lib$(TARGET).a : ${OBJ}
# 	${AR} -r $@ ${OBJ}

lib$(TARGET).so : ${OBJ}
	$(CXX) $(CXXFLAGS) ${LDFLAGS} -o lib$(TARGET).so ${OBJ} $(SRC) $(LIBS)   
	# -o $@ ${OBJECTS_MCMC}

# $(TARGET): ${OBJECTS_MCMC}
# 	$(CXX) $(CXXFLAGS) -o $(TARGET) ${OBJECTS_MCMC} $(LIBS)

# cleanup and install
.PHONY: clean
clean:
	@rm -rf *.so ./*.gcov ./*.gcno ./*.gcda ./*.dSYM ./*.test ./*.o ./*.a


# from fabber toymodel
# #
# # Build
# #

# all:	${XFILES} libfabber_models_${MODELS_NAME}.a

# # models in a library
# libfabber_models_${MODELS_NAME}.a : ${OBJS}
# 	${AR} -r $@ ${OBJS}

# # fabber built from the FSL fabbercore library including the models specifieid in this project
# fabber_${MODELS_NAME} : fabber_main.o ${OBJS}
# 	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${OBJS} -lfabbercore -lfabberexec ${LIBS}

# # DO NOT DELETE