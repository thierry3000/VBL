

find_package(OpenMP)
if(OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
find_package(Boost REQUIRED COMPONENTS "random")

add_executable(test_omp_slurm test_omp_slurm.cpp)
target_link_libraries(test_omp_slurm boost_random)
