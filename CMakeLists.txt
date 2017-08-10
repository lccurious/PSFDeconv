project(psf)

cmake_minimum_required(VERSION 2.9)

set(CMAKE_CXX_STANDARD 11)

include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin )
set(SOURCE_FILES src/main.cpp
		include/comparemat.h
		src/comparemat.cpp
		include/cvutils.h
		src/cvutils.cpp
		include/deconvpsf.h
		src/deconvpsf.cpp
        include/psfgen.h
        src/psfgen.cpp
		include/vtkutils.h
		src/vtkutils.cpp
		include/config.h)

# BOOST LIB CONFIG
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package( Boost REQUIRED COMPONENTS thread program_options)


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
	add_executable( psf ${SOURCE_FILES} )
	target_link_libraries( psf ${Boost_LIBRARIES} matio fftw3 tiff)
#	add_executable( psf2Ddemo src/psfg.cpp
#			src/psf_gen2.cpp
#			src/cvTools.cpp)
#	target_link_libraries( psf2Ddemo ${Boost_LIBRARIES} )
#
#	add_executable( psfSpeed src/psf_accelerate.cpp
#			src/psf_gen2.cpp )
#	target_link_libraries( psfSpeed ${Boost_LIBRARIES} )
#
#	add_executable( deconvLucyDemo src/demo_deconv_LR.cpp )
endif()

if(OpenCV_LIBS)
#	target_link_libraries( psf2Ddemo ${OpenCV_LIBS})
#	target_link_libraries( deconvLucyDemo ${OpenCV_LIBS} )
    target_link_libraries( psf ${OpenCV_LIBS} )
endif()

if(VTK_LIBRARIES)
    target_link_libraries( psf ${VTK_LIBRARIES} )
#	target_link_libraries( psf2Ddemo ${VTK_LIBRARIES})
endif()
