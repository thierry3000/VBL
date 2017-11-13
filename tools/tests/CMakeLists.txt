add_executable(check_random check_random.cpp)
find_package(TBB REQUIRED)
message("TBB_LIBRARIES_DIR: ${TBB_LIBRARIES_DIR}")
message("TBB_LIBRARY: ${TBB_LIBRARY}")
#TF: somehow the cmake is broken here, I overwrite with environmental variabels
message(STATUS "TBB_LIBRARIES          	= $ENV{TBB_LIBRARY_RELEASE}")
set(TBB_LIBRARIES_DIR $ENV{TBB_LIBRARY_RELEASE})
message(STATUS "TBB_INCLUDE_DIR 	      = $ENV{TBB_INCLUDE_DIR}")
set(TBB_INCLUDE_DIRS $ENV{TBB_INCLUDE_DIR})
add_executable(test_delauny test_delauny.cpp)
target_link_libraries(test_delauny gmp CGAL tbb)

add_executable(test_delauny_parallel test_delauny_parallel.cpp)
target_link_libraries(test_delauny_parallel gmp CGAL tbb)
