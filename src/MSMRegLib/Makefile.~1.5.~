# $id: Makefile,v 1.1.1.1 2013/10/08 14:06:18 emmar Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = MeshReg
USRCXXFLAGS=-g 

USRINCFLAGS =-I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_ZLIB} -I${INC_CPROB} -I${INC_PROB}  -I${INC_GD}  -I${INC_GDC} -I${INC_PNG} -I${INC_BOOST} -I${HOME}/fsldev/include -I${HOME}/fsldev/extras/include
USRLDFLAGS =-L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_ZLIB}  -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_GD}   -L. -L${HOME}/fsldev/lib  -L${HOME}/extras/fsldev/lib
LIBS= -lmsm -lfastpd -ldiscrete -lnewmeshclass -lfslsurface -lgiftiio -lmm  -lexpat  -lfslvtkio -lnewimage  -lmiscmaths -lmiscplot -lmiscpic -lfslio -lniftiio -lnewmat -lutils  -lbint -lnewran -lprob -lznz -lgdc -lgd -lpng  -lz
MSMLIBOBJS= ContinuosOpt.o meshmodify.o meshreg.o

ifdef HOCR
	USRCXXFLAGS:= -std=c++11  -DHAS_HOCR
endif

ifdef GROUP
	MYPROGRAMOBJ:= groupmeshreg.o   	
endif

ifdef TBB
	USRCXXFLAGS:= -DHAS_TBB
	LIBS= -lmsm -lfastpd -ldiscrete -ltbb  -lANN -lnewmeshclass -lfslsurface -lgiftiio -lmm  -lexpat -lfirst_lib -lshapeModel -lfslvtkio -lmeshclass -lnewimage  -lmiscmaths -lmiscplot -lmiscpic -lfslio -lniftiio -lnewmat -lutils  -lbint -lnewran -lprob -lznz -lgdc -lgd -lpng  -lz
endif 




all: libmeshreg.a 

libmeshreg.a: ${MSMLIBOBJS}
	${AR} -r libmeshreg.a ${MSMLIBOBJS}


