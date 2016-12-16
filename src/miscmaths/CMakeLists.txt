SET(MISCMATHS_SOURCES
  miscmaths.cc
  optimise.cc
  miscprob.cc
  kernel.cc
  histogram.cc
  base2z.cc
  t2z.cc
  f2z.cc
  minimize.cc
  cspline.cc
  sparse_matrix.cc
  sparsefn.cc
  rungekutta.cc)

INCLUDE_DIRECTORIES( ${FSL_SOURCE_DIR}/libprob )

ADD_LIBRARY(miscmaths STATIC ${MISCMATHS_SOURCES})

IF(FSL_BUILD)
  INSTALL_FILES(/include/miscmaths .h ${MISCMATHS_SOURCES})
  INSTALL_TARGETS(/lib miscmaths)
ENDIF(FSL_BUILD)
