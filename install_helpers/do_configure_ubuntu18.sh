cmake \
    -DCMAKE_INSTALL_PREFIX=/localdisk/vbl_install \
    -DCMAKE_BUILD_TYPE=Debug \
    -DADDITIONAL_INCLUDE_DIRS="" \
    -DADDITIONAL_LIBRARY_DIRS="" \
    -DHDF5_INCLUDE_DIRS="/usr/include/hdf5/openmpi" \
    -DCMAKE_CXX_COMPILER="/usr/bin/mpic++" \
    -DCGAL_DIR="/usr/lib/x86_64-linux-gnu/cmake/CGAL" \
    -DKREBSLIB="/localdisk/tc_install/lib" \
    $1
