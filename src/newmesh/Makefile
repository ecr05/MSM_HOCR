# $Id: Makefile,v 1.14 2016/05/12 12:29:42 cowboy Exp $
include ${FSLCONFDIR}/default.mk

# changed projectname from newmeshclass for compatability
PROJNAME = newmesh

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB} -I${FSLDIR}/include/niftiio  -I${INC_BOOST} 
USRLDFLAGS = -L${LIB_NEWMAT}  -L${LIB_ZLIB} -L./

LIBS=-lfslsurface -lgiftiio  -lexpat -lfirst_lib -lshapeModel -lfslvtkio -lmeshclass -lnewimage  -lmiscmaths -lfslio -lniftiio -lnewmat -lutils  -lbint -lnewran -lprob -lznz -lz

LIBFILES=libnewmeshclass.a
XFILES=examples

MESHOBJS= giftiInterface.o point.o newmesh.o Relations.o resampler.o meshfns.o featurespace.o  

all: libnewmeshclass.a 

libnewmeshclass.a: ${MESHOBJS}
	${AR} -r libnewmeshclass.a ${MESHOBJS} 



