rm -fr build
mkdir build
cd build
BZ=-DWITH_COMPRESSION=bzip2
UL=/usr/local
PPATH="$UL"
#NONC4="-DENABLE_NETCDF_4=OFF"
HDF5="-DHDF5_LIB=${UL}/lib/libhdf5.so -DHDF5_HL_LIB=${UL}/lib/libhdf5_hl.so -DHDF5_INCLUDE_DIR=${UL}/include"
cmake -DCMAKE_INSTALL_PREFIX=${UL} -DCMAKE_PREFIX_PATH="$PPATH" ${HDF5} ${NONC4} ${BZ} ..
cmake --build .
make test
