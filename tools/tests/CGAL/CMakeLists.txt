
find_package(TBB REQUIRED)
#message("TBB_LIBRARIES_DIR: ${TBB_LIBRARIES_DIR}")
#message("TBB_LIBRARY: ${TBB_LIBRARY}")
# #T.F.: somehow the cmake is broken here, I overwrite with environmental variabels
# message(STATUS "TBB_LIBRARIES          	= $ENV{TBB_LIBRARY_RELEASE}")
# set(TBB_LIBRARIES_DIR $ENV{TBB_LIBRARY_RELEASE})
# message(STATUS "TBB_INCLUDE_DIR 	      = $ENV{TBB_INCLUDE_DIR}")
# set(TBB_INCLUDE_DIRS $ENV{TBB_INCLUDE_DIR})



# add_executable(test_delauny_parallel test_delauny_parallel.cpp)
# target_link_libraries(test_delauny_parallel gmp CGAL)
add_executable(test_delauny test_delauny.cpp)
target_link_libraries(test_delauny "${myLibs}")
add_executable(test_delauny_parallel test_delauny_parallel.cpp)
target_link_libraries(test_delauny_parallel "${myLibs}")


add_executable(test_regular_parallel test_regular_triangulation_parallel.cpp)
target_link_libraries(test_regular_parallel ${myLibs})
add_executable(info_insert_with_pair_iterator_parallel info_insert_with_pair_iterator_parallel.cpp)
target_link_libraries(info_insert_with_pair_iterator_parallel ${myLibs})

add_executable(info_insert_with_pair_iterator info_insert_with_pair_iterator.cpp)
target_link_libraries(info_insert_with_pair_iterator ${myLibs})
#note: tbb is not used in this example but required since CGAL is linked against it
#target_link_libraries(info_insert_with_pair_iterator ${TBB_LIBRARIES})
add_executable(parallel_insertion_in_delaunay_3 parallel_insertion_in_delaunay_3.cpp)
target_link_libraries(parallel_insertion_in_delaunay_3 ${myLibs})

add_executable(parallel_insertion_in_delaunay_3b parallel_insertion_in_delaunay_3b.cpp)
target_link_libraries(parallel_insertion_in_delaunay_3b ${myLibs})

add_executable(max_neighbor max_neighbor.cpp)
target_link_libraries(max_neighbor ${myLibs})

add_executable(linear_interpolation linear_interpolation.cpp)
target_link_libraries(linear_interpolation ${myLibs})
