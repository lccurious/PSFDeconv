project(psf)

cmake_minimum_required(VERSION 2.9)

set(CMAKE_CXX_STANDARD 11)

# BOOST LIB CONFIG
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package( Boost REQUIRED COMPONENTS thread )

#OPENCV CONFIG
find_package( OpenCV REQUIRED )
if(OpenCV_LIBS)
	include_directories( ${OpenCV_INCLUDE_DIRS} )
endif()

#VTK CONFIG
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
	add_executable( bin/psf_2Ddemo src/psfg.cpp )
	target_link_libraries( bin/psf_2Ddemo ${Boost_LIBRARIES} )
	add_executable( bin/psf_speed src/psf_accelerate.cpp )
	target_link_libraries( bin/psf_speed ${Boost_LIBRARIES} )
endif()

if(OpenCV_LIBS)
	target_link_libraries( bin/psf_2Ddemo ${OpenCV_LIBS})
endif()

if(VTK_LIBRARIES)
	target_link_libraries( bin/psf_2Ddemo ${VTK_LIBRARIES})
endif()
