

add_library(vbl SHARED Environment.cpp EnvironmentalSignals.cpp CellType.cpp CellsSystem-A.cpp CellsSystem-B.cpp CellsSystem-C.cpp CellsSystem-D-2.cpp CellsSystem-E.cpp CellsSystem-F.cpp Utilities.cpp BloodVessel.cpp CellsSystem.cpp)
add_executable(vbl-v3.2 main.cpp )

message("linking to: ${myLibs}")
#target_link_libraries(vbl "${myLibs}")
message("blubl")
set_target_properties(vbl PROPERTIES INSTALL_RPATH_USE_LINK_PATH FALSE)
target_link_libraries(vbl boost_thread CGAL gmp gmpxx)
target_link_libraries(vbl-v3.2 vbl)

install(DIRECTORY . DESTINATION include/vbl FILES_MATCHING PATTERN "*.h")
install(TARGETS vbl LIBRARY DESTINATION lib)
install(TARGETS vbl-v3.2 DESTINATION bin)
