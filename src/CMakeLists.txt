

add_library(vbl SHARED Environment.cpp EnvironmentalSignals.cpp CellType.cpp CellsSystem-A.cpp CellsSystem-B.cpp CellsSystem-C.cpp CellsSystem-D-2.cpp CellsSystem-E.cpp CellsSystem-F.cpp Utilities.cpp BloodVessel.cpp CellsSystem.cpp)
set_target_properties(vbl PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(vbl PROPERTIES SOVERSION 1)

set_target_properties(vbl PROPERTIES PUBLIC_HEADER ../include/vbl.h)

message("linking to: ${myLibs}")
#target_link_libraries(vbl "${myLibs}")

target_link_libraries(vbl "${myLibs}")
set_target_properties(vbl PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(vbl-v3.2 main.cpp )
target_include_directories(vbl-v3.2 PRIVATE ../include)
target_link_libraries(vbl-v3.2 vbl)
set_target_properties(vbl-v3.2 PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

### NOTE THIS IS A HACK AND SHOULD BE DONE PROPERLY IN FUTURE!!!
message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};/localdisk/thierry/local/lib/CGAL;/localdisk/thierry/local/lib64/CGAL")
message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
if( ACTIVATE_CONCURRENCY )
  find_package( TBB REQUIRED )
  message("TBB_LIBRARIES: ${TBB_LIBRARIES}") # variable used by CGAL
  if( TBB_FOUND )
    message("tbb found!")
    include(${TBB_USE_FILE})
    target_link_libraries(vbl ${TBB_LIBRARIES})
    target_link_libraries(vbl-v3.2 ${TBB_LIBRARIES})
  endif()
endif()

install(TARGETS vbl-v3.2 DESTINATION bin)

message("header goes to: ${CMAKE_INSTALL_INCLUDEDIR}")
install(TARGETS vbl
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include)

install(TARGETS vbl LIBRARY DESTINATION lib)
