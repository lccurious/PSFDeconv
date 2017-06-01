project(cell)
cmake_minimum_required(VERSION 2.9)

find_package(OpenCV REQUIRED)

message(STATUS "OpenCV library status:")
message(STATUS "    version: ${OpenCV_VERSION}")
message(STATUS "    libraries: ${OpenCV_LIBS}")
message(STATUS "    include path: ${OpenCV_INCLUDE_DIRS}")

include_directories( ${OpenCV_INCLUDE_DIRS} )

add_executable( bin/cell src/conv_main.cpp)

target_link_libraries( bin/cell ${OpenCV_LIBS} )
