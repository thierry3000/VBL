add_library(milotti_mts SHARED Environment.cpp EnvironmentalSignals.cpp CellType.cpp CellsSystem.cpp CellsSystem-A.cpp CellsSystem-B.cpp CellsSystem-C.cpp CellsSystem-D-2.cpp CellsSystem-E.cpp CellsSystem-F.cpp Utilities.cpp BloodVessel.cpp)
add_executable(Sim3D-v3.2 main.cpp )
#add_executable(Sim3D-v3.2 main.cpp Environment.cpp EnvironmentalSignals.cpp CellType.cpp CellsSystem.cpp CellsSystem-D-2.cpp Utilities.cpp)
target_link_libraries(Sim3D-v3.2 "milotti_mts;${myLibs}")

install(DIRECTORY . DESTINATION milotti_mts/include FILES_MATCHING PATTERN "*.h")
install(TARGETS milotti_mts LIBRARY DESTINATION milotti_mts/lib)
install(TARGETS Sim3D-v3.2 DESTINATION milotti_mts/bin)
