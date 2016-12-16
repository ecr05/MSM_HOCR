# A Makefile for the miscmaths class

include ${FSLCONFDIR}/default.mk

PROJNAME = miscmaths

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_BOOST} -I${INC_PROB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB}

OBJS = miscmaths.o optimise.o miscprob.o kernel.o histogram.o base2z.o t2z.o f2z.o minimize.o cspline.o sparse_matrix.o sparsefn.o rungekutta.o nonlin.o bfmatrix.o Simplex.o

LIBS = -lutils -lnewmat -lprob -lm

# The target "all" should ALWAYS be provided
#   typically it will just be another target name

all: libmiscmaths.a

quick:${OBJS} quick.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} quick.o ${LIBS}

libmiscmaths.a: ${OBJS}
	${AR} -r libmiscmaths.a ${OBJS}
