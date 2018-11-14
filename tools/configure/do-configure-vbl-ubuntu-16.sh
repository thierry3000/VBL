#note old gentoo sytem: -has no CImg
#                       - trilinos is not relative to eigen, so need eigen explicitly

cmake \
    -D ACTIVATE_CONCURRENCY=ON \
    -D ADDITIONAL_INCLUDE_DIRS="/localdisk/tc_install/include" \
    -D ADDITIONAL_LIBRARY_DIRS="/localdisk/tc_install/lib" \
    -D KREBSLIB="/localdisk/tc_install/lib/libkrebs_.so" \
    -D CMAKE_INSTALL_PREFIX="/usr/local" \
    -D CMAKE_BUILD_TYPE=Release \
    -D CGAL_DIR="/usr/lib/x86_64-linux-gnu/cmake/CGAL" \
    -D HDF5_INCLUDE_DIRS=/usr/lib/x86_64-linux-gnu/hdf5/serial/include \
    $1
