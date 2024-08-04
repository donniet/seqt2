cmake_minimum_required(VERSION 3.5.0)
project(seqt2 VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)

find_package(Boost 1.74.0 REQUIRED COMPONENTS log)
find_package(OpenCL REQUIRED)

include_directories(include)

add_library(seqt src/seqt.cpp)
target_compile_definitions(seqt PRIVATE BOOST_COMPUTE_DEBUG_KERNEL_COMPILATION)
target_link_libraries(seqt ${BOOST_LIBRARIES} ${OpenCL_LIBRARIES})

add_executable(main main.cpp)
target_link_libraries(main seqt)

