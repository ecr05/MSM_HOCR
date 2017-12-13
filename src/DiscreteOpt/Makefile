
include ${FSLCONFDIR}/default.mk

PROJNAME=DiscreteOpt
HOCR=True
USRCXXFLAGS=-g -DNDEBUG  
USRINCFLAGS =-I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_ZLIB} -I${INC_CPROB} -I${INC_PROB} -I${INC_BOOST} -I${FSLDEVDIR}/include -I${FSLDEVDIR}/extras/include
USRLDFLAGS =-L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_ZLIB} -L${LIB_CPROB} -L${LIB_PROB} -L${FSLDEVDIR}/lib #${TBBLIBPATH}

LIBS=-lnewmeshclass -lfslsurface -lgiftiio  -lexpat -lfslvtkio -lnewimage  -lmiscmaths -lfslio -lniftiio -lnewmat -lutils  -lbint -lnewran -lprob -lznz -lz
MYPROGRAMOBJ=similarities.o DiscreteCostFunction.o  DiscreteModel.o 
MYPROGRAM=DiscreteOpt
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
	LIBS=-ltbb -lnewmeshclass -lfslsurface -lgiftiio  -lexpat -lfslvtkio -lnewimage  -lmiscmaths -lfslio -lniftiio -lnewmat -lutils  -lbint -lnewran -lprob -lznz -lz
endif
   


#-----------------------------------------------------------------		

MYINCLUDES=

MYLIBRARIES=

#------------------------------------------------------------------------------

all: ${MYlib}.a

${MYlib}.a: ${MYPROGRAMOBJ}
	${AR} -r libdiscrete.a ${MYPROGRAMOBJ}

