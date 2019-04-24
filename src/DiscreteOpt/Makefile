
include ${FSLCONFDIR}/default.mk

PROJNAME=DiscreteOpt
USRCXXFLAGS=-g -DNDEBUG
USRINCFLAGS =-I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_ZLIB} -I${INC_CPROB} -I${INC_PROB} -I${INC_BOOST} -I${FSLDEVDIR}/include -I${FSLDEVDIR}/extras/include
USRLDFLAGS =-L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_ZLIB} -L${LIB_CPROB} -L${LIB_PROB} -L${FSLDEVDIR}/lib #${TBBLIBPATH}

MYPROGRAMOBJ=similarities.o DiscreteCostFunction.o  DiscreteModel.o
Mylib=libdiscrete

ifdef HOCR
	USRCXXFLAGS+= -std=c++11  -DHAS_HOCR  -DNDEBUG
	MSMLIBOBJS+= Fusion.o
endif

ifdef GROUP
	MYPROGRAMOBJ+= DiscreteGroupCostFunction.o  DiscreteGroupModel.o
endif

ifdef TBB
	USRCXXFLAGS+= -DHAS_TBB
endif

all: ${MYlib}.a

${MYlib}.a: ${MYPROGRAMOBJ}
	${AR} -r libdiscrete.a ${MYPROGRAMOBJ}
