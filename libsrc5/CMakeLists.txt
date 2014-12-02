
# Process these files with m4.
SET(m4_SOURCES nc5attr nc5ncx nc5putget)
foreach (f ${m4_SOURCES})      
  GEN_m4(${f})
endforeach(f)

SET(libsrc5_SOURCES nc5v1hpg.c nc5putget.c nc5attr.c nc5dispatch.c nc5internal.c nc5var.c nc5dim.c nc5ncx.c)

add_library(netcdf5 OBJECT ${libsrc5_SOURCES}) 

IF (USE_PNETCDF)
  SET(ARGS_MANPAGE ${ARGS_MANPAGE} -DUSE_PNETCDF=TRUE)
ENDIF(USE_PNETCDF)

IF (BUILD_PARALLEL)
  SET(ARGS_MANPAGE ${ARGS_MANPAGE} -DPARALLEL_IO=TRUE)
ENDIF(BUILD_PARALLEL)

ADD_EXTRA_DIST(${libsrc5_SOURCES} ${m4_SOURCES} nc5ncx.h CMakeLists.txt Makefile.am)
