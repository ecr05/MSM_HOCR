include ${FSLCONFDIR}/default.mk

PROJNAME= MSM
USRCXXFLAGS=-g
HOCR=True

USRINCFLAGS=-I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_ZLIB} -I${INC_CPROB} -I${INC_PROB}  -I${INC_GD}  -I${INC_GDC} -I${INC_PNG} -I${INC_BOOST} -I./include -I${FSLDEVDIR}/include -I${FSLDEVDIR}/extras/include -I${FSLDEVDIR}/extras/ELC1.04
USRLDFLAGS=-L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_ZLIB}  -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_GD}   -L. -L${FSLDEVDIR}/lib  -L${FSLDEVDIR}/extras/lib
LIBS= -lmeshreg -lfastpd -ldiscrete -lnewmeshclass -lfslsurface -lgiftiio -lmm -lexpat -lfslvtkio -lnewimage -lmiscmaths -lmiscplot -lmiscpic -lNewNifti -lutils  -lbint -lnewran -lprob -lznz -lgdc -lgd -lpng  -lz

ifdef HOCR
	USRCXXFLAGS+= -std=c++11  -DHAS_HOCR
endif

ifdef TBB
	USRCXXFLAGS+= -DHAS_TBB
	LIBS= -lmeshreg -lfastpd -ldiscrete -ltbb -lnewmeshclass -lfslsurface -lgiftiio -lmm -lexpat -lfslvtkio -lnewimage -lmiscmaths -lmiscplot -lmiscpic -lNewNifti -lutils  -lbint -lnewran -lprob -lznz -lgdc -lgd -lpng  -lz

endif

ifdef GROUP
	all:= msm_group msm_metricmath msmresamplemesh  estimate_strains estimate_curvature msm_metric_average msm_makeicos intensity_normalise minimise_anat-to-sphere_distortions:
endif

all:msm  msmapplywarp  estimate_metric_distortion surfconvert msmresample msm_metric_sim

libmsm.a: ${MSMLIBOBJS}
	${AR} -r libmsm.a ${MSMLIBOBJS}


msm: msm.o
	${CXX} ${LDFLAGS} -o $@ msm.o ${LIBS}
msm_group: msm_group.o
	${CXX} ${LDFLAGS} -o $@ msm_group.o ${LIBS}

msm_makeicos: msm_makeicos.o
	${CXX} ${LDFLAGS} -o $@ msm_makeicos.o ${LIBS}
intensity_normalise: intensity_normalise.o
	${CXX} ${LDFLAGS} -o $@ intensity_normalise.o ${LIBS}
msm_metricmath: msm_metricmath.o
	${CXX} ${LDFLAGS} -o $@ msm_metricmath.o ${LIBS}
msm_surface_average: msm_surface_average.o
	${CXX} ${LDFLAGS} -o $@ msm_surface_average.o ${LIBS}
msm_metric_sim: msm_metric_sim.o
	${CXX} ${LDFLAGS} -o $@ msm_metric_sim.o ${LIBS}

estimate_strains: estimate_strains.o
	${CXX} ${LDFLAGS} -o $@ estimate_strains.o ${LIBS}
estimate_curvature: estimate_curvature.o
	${CXX} ${LDFLAGS} -o $@ estimate_curvature.o ${LIBS}

groupsegmentation: groupsegmentation.o
	${CXX} ${LDFLAGS} -o $@ groupsegmentation.o ${LIBS}

msmapplywarp: msmapplywarp.o
	${CXX}  ${LDFLAGS} -o $@  msmapplywarp.o ${LIBS}

msmresample: msmresample.o
	${CXX}  ${LDFLAGS} -o $@  msmresample.o ${LIBS}

msmresamplemesh: msmresamplemesh.o
	${CXX} ${LDFLAGS} -o $@ msmresamplemesh.o ${LIBS}

estimate_metric_distortion: estimate_metric_distortion.o
	${CXX} ${LDFLAGS} -o $@ estimate_metric_distortion.o ${LIBS}

minimise_anat-to-sphere_distortions: minimise_anat-to-sphere_distortions.o
	${CXX} ${LDFLAGS} -o $@ minimise_anat-to-sphere_distortions.o ${LIBS}

msm_metric_average: msm_metric_average.o
	${CXX}  ${LDFLAGS} -o $@  msm_metric_average.o ${LIBS}

surfconvert: surfconvert.o
	${CXX}  ${LDFLAGS} -o $@  surfconvert.o ${LIBS}
