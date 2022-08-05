include ${FSLCONFDIR}/default.mk

PROJNAME = MeshReg
USRCXXFLAGS=-g
HOCR=True

USRINCFLAGS =-I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_ZLIB} -I${INC_CPROB} -I${INC_PROB}  -I${INC_GD}  -I${INC_GDC} -I${INC_PNG} -I${INC_BOOST} -I${FSLDEVDIR}/include -I${FSLDEVDIR}/extras/include -I${FSLDEVDIR}/extras/ELC1.04
USRLDFLAGS =-L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_ZLIB}  -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_GD}   -L. -L${FSLDEVDIR}/lib  -L${FSLDEVDIR}/extras/lib
MSMLIBOBJS= ContinuosOpt.o meshmodify.o meshreg.o

ifdef HOCR
	USRCXXFLAGS+= -std=c++11  -DHAS_HOCR  -DNDEBUG
	MSMLIBOBJS+= Fusion.o
endif

ifdef GROUP
	MSMLIBOBJS+= groupmeshreg.o
endif

ifdef TBB
	USRCXXFLAGS+= -DHAS_TBB
endif




all: libmeshreg.a

libmeshreg.a: ${MSMLIBOBJS}
	${AR} -r libmeshreg.a ${MSMLIBOBJS}
