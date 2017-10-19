

add_library(vbl SHARED Environment.cpp EnvironmentalSignals.cpp CellType.cpp CellsSystem-A.cpp CellsSystem-B.cpp CellsSystem-C.cpp CellsSystem-D-2.cpp CellsSystem-E.cpp CellsSystem-F.cpp Utilities.cpp BloodVessel.cpp CellsSystem.cpp)
set_target_properties(vbl PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(vbl PROPERTIES SOVERSION 1)

set_target_properties(vbl PROPERTIES PUBLIC_HEADER ../include/vbl.h)

message("linking to: ${myLibs}")
#target_link_libraries(vbl "${myLibs}")

target_link_libraries(vbl boost_thread CGAL gmp gmpxx)


# add_executable(vbl-v3.2 main.cpp )
# target_include_directories(vbl-v3.2 PRIVATE ../include)
# target_link_libraries(vbl-v3.2 vbl)
# install(TARGETS vbl-v3.2 DESTINATION bin)

message("header goes to: ${CMAKE_INSTALL_INCLUDEDIR}")
install(TARGETS vbl
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include)

install(TARGETS vbl LIBRARY DESTINATION lib)
