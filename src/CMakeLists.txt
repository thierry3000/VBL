add_library(vbl SHARED Environment.cpp EnvironmentalSignals.cpp CellType.cpp CellsSystem.cpp CellsSystem-A.cpp CellsSystem-B.cpp CellsSystem-C.cpp CellsSystem-D-2.cpp CellsSystem-E.cpp CellsSystem-F.cpp Utilities.cpp BloodVessel.cpp)
add_executable(vbl-v3.2 main.cpp )
#add_executable(Sim3D-v3.2 main.cpp Environment.cpp EnvironmentalSignals.cpp CellType.cpp CellsSystem.cpp CellsSystem-D-2.cpp Utilities.cpp)
target_link_libraries(vbl-v3.2 "milotti_mts;${myLibs}")

install(DIRECTORY . DESTINATION include/vbl FILES_MATCHING PATTERN "*.h")
install(TARGETS vbl LIBRARY DESTINATION lib/vbl)
install(TARGETS vbl-v3.2 DESTINATION bin/vbl)
