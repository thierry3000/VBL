add_executable(check_random check_random.cpp)
target_link_libraries(check_random boost_random)
# find_package(TBB REQUIRED)
# message("TBB_LIBRARIES_DIR: ${TBB_LIBRARIES_DIR}")
# message("TBB_LIBRARY: ${TBB_LIBRARY}")
# #TF: somehow the cmake is broken here, I overwrite with environmental variabels
# message(STATUS "TBB_LIBRARIES          	= $ENV{TBB_LIBRARY_RELEASE}")
# set(TBB_LIBRARIES_DIR $ENV{TBB_LIBRARY_RELEASE})
# message(STATUS "TBB_INCLUDE_DIR 	      = $ENV{TBB_INCLUDE_DIR}")
# set(TBB_INCLUDE_DIRS $ENV{TBB_INCLUDE_DIR})
add_executable(test_delauny test_delauny.cpp)
target_link_libraries(test_delauny ${myLibs})

# add_executable(test_delauny_parallel test_delauny_parallel.cpp)
# target_link_libraries(test_delauny_parallel gmp CGAL)

add_executable(test_regular_parallel test_regular_triangulation_parallel.cpp)
target_link_libraries(test_regular_parallel ${myLibs})

#receive by email from Edoardo Milotti 14.11.2017
add_executable(test_thread test_thread.cpp)
target_link_libraries(test_thread ${myLibs})

message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};/localdisk/thierry/local/lib/CGAL;")
message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

if( ACTIVATE_CONCURRENCY )
  find_package( TBB REQUIRED )
  message("TBB_LIBRARIES: ${TBB_LIBRARIES}") # variable used by CGAL
  if( TBB_FOUND )
    message("tbb found!")
    message("TBB_USE_FILE: ${TBB_USE_FILE}")
    include(${TBB_USE_FILE})
    add_executable(test_delauny_parallel test_delauny_parallel.cpp)
    #target_link_libraries(test_delauny_parallel "${CGAL_LIBRARIES};${TBB_LIBRARIES}")
    target_link_libraries(test_delauny_parallel ${myLibs})
    target_link_libraries(test_delauny_parallel ${TBB_LIBRARIES})
    target_link_libraries(test_delauny ${TBB_LIBRARIES})
    target_link_libraries(test_thread ${TBB_LIBRARIES})
    target_link_libraries(test_regular_parallel ${TBB_LIBRARIES})
    add_executable(info_insert_with_pair_iterator_parallel info_insert_with_pair_iterator_parallel.cpp)
    target_link_libraries(info_insert_with_pair_iterator_parallel ${myLibs})
    target_link_libraries(info_insert_with_pair_iterator_parallel ${TBB_LIBRARIES})
    add_executable(info_insert_with_pair_iterator info_insert_with_pair_iterator.cpp)
    target_link_libraries(info_insert_with_pair_iterator ${myLibs})
    #note: tbb is not used in this example but required since CGAL is linked against it
    target_link_libraries(info_insert_with_pair_iterator ${TBB_LIBRARIES})
  endif()
endif()




