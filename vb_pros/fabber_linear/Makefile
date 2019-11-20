include ${FSLCONFDIR}/default.mk

MODELS_NAME=linear1

PROJNAME = fabber_${MODELS_NAME}

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST} -I..
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L/lib64


FSLVERSION= $(shell cat ${FSLDIR}/etc/fslversion | head -c 1)
ifeq ($(FSLVERSION), 5) 
  NIFTILIB = -lfslio -lniftiio 
else 
  NIFTILIB = -lNewNifti
endif

LIBS = -lutils -lnewimage -lmiscmaths -lprob -lnewmat ${NIFTILIB} -lznz -lz -ldl
TESTLIBS = -lgtest -lpthread

XFILES = fabber_${MODELS_NAME}

# Forward models
OBJS =  fwdmodel_${MODELS_NAME}.o

# For debugging:
OPTFLAGS = -ggdb
#OPTFLAGS =

#
# Build
#

all:	${XFILES} libfabber_models_${MODELS_NAME}.a

# models in a library
libfabber_models_${MODELS_NAME}.a : ${OBJS}
	${AR} -r $@ ${OBJS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_${MODELS_NAME} : fabber_main.o ${OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${OBJS} -lfabbercore -lfabberexec ${LIBS}

# DO NOT DELETE
