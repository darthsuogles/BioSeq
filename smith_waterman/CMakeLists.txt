cmake_minimum_required(VERSION 2.6.2)

project (alignment)
set (alignment_VERSION_MAJOR 0)
set (alignment_VERSION_MINOR 1)

find_package(CUDA REQUIRED)
find_package(doxygen REQUIRED)


if( CUDA_CUDA_LIBRARY)
  message(STATUS "CUDA found on this system")
  
  # this block sets up the shared library called proj1
  # set (CUDA_NVCC_FLAGS --shared  --compiler-options;'-fPIC' )
  # include_directories(${CUDA_INCLUDE_DIRS}) 
  # cuda_add_library(prog1 SHARED cu_lib.cu)
  # target_link_libraries(prog1 ${CUDA_LIBRARIES})
  
  # this compiles and links the binary, which is called align
  include_directories(${CMAKE_SOURCE_DIR}) 
  add_executable(align sw_align_cpu.cpp)
  set (libs_to_link ${CUDA_LIBRARIES})
  target_link_libraries(align ${libs_to_link})  

else( CUDA_CUDA_LIBRARY)
  message(STATUS "oops CUDA not found on this system")
endif(CUDA_CUDA_LIBRARY)

# This set all the optimization flags
set ( CMAKE_BUILD_TYPE Release )
