

#add_library(vbl SHARED Environment.cpp EnvironmentalSignals.cpp CellType.cpp CellsSystem-A.cpp CellsSystem-B.cpp CellsSystem-C.cpp CellsSystem-D.cpp CellsSystem-E.cpp CePROJECT_SOURCE_DIRllsSystem-F.cpp Utilities.cpp BloodVessel.cpp CellsSystem.cpp)
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/include")
add_library(vbl SHARED Environment.cpp EnvironmentalSignals.cpp CellType.cpp Utilities.cpp BloodVessel.cpp CellsSystem.cpp)
#set_target_properties(vbl PROPERTIES VERSION ${PROJECT_VERSION})
#set_target_properties(vbl PROPERTIES SOVERSION 1)

#set_target_properties(vbl PROPERTIES PUBLIC_HEADER ../include/vbl.h)
if(USE_TUMORCODE)
  message("ADDITIONAL_LIBRARY_DIRS: ${ADDITIONAL_LIBRARY_DIRS}")
  if(CMAKE_BUILD_TYPE MATCHES Debug)
    find_library(KREBSLIB NAMES krebs_d PATHS ${ADDITIONAL_LIBRARY_DIRS} )
  endif()
  if(CMAKE_BUILD_TYPE MATCHES Release)
    find_library(KREBSLIB NAMES krebs_ PATHS ${ADDITIONAL_LIBRARY_DIRS} )
  endif()
  message("KREBSLIB: ${KREBSLIB}")
  find_package(HDF5 REQUIRED COMPONENTS CXX)
  include_directories(${HDF5_INCLUDE_DIRS})
  set(myLibs "${myLibs};${KREBSLIB};${HDF5_LIBRARY_DIRS}")
endif()

message("linking to: ${myLibs}")
#target_link_libraries(vbl "${myLibs}")

target_link_libraries(vbl "${myLibs}")
target_link_libraries(vbl gmp)
#set_target_properties(vbl PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(vbl_cells main.cpp )
target_include_directories(vbl_cells PRIVATE ../include)
target_link_libraries(vbl_cells vbl)
#set_target_properties(vbl_cells PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS vbl_cells DESTINATION bin)

message("header goes to: ${CMAKE_INSTALL_INCLUDEDIR}")
install(TARGETS vbl
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include)

install(TARGETS vbl LIBRARY DESTINATION lib)
