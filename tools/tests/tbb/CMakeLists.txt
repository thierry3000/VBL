

if( ACTIVATE_CONCURRENCY )
  find_package( TBB REQUIRED )
  message("TBB_LIBRARIES: ${TBB_LIBRARIES}") # variable used by CGAL
  if( TBB_FOUND )
    message("tbb found!")
    message("TBB_USE_FILE: ${TBB_USE_FILE}")
    include(${TBB_USE_FILE})
    add_executable(test_tbb test_tbb_for.cpp)
    target_link_libraries(test_tbb ${myLibs})
    target_link_libraries(test_tbb ${TBB_LIBRARIES})
    add_executable(parallel_do parallel_do.cpp)
    target_link_libraries(parallel_do ${myLibs})
    target_link_libraries(parallel_do ${TBB_LIBRARIES})
  endif()
endif()
