cmake_minimum_required(VERSION 2.6)

project(libpointmesher)

#include(UseDoxygen)

# boost
find_package(Boost REQUIRED)

# eigen 2
find_path(EIGEN2_INCLUDE_DIR Eigen/Core
	/usr/local/include/eigen2
	/usr/include/eigen2
)

# nabo
find_path(NABO_INCLUDE_DIR nabo/nabo.h
	/usr/local/include
	/usr/include
)
find_library(NABO_LIBRARY libnabo
	/usr/local/lib
	/usr/lib
)

# look for OpenCL library
if (WIN32)
	find_library(OPENCL_LIBRARIES opencl64)
	if (!OPENCL_LIBRARIES)
		find_library(OPENCL_LIBRARIES opencl32)
	endif (!OPENCL_LIBRARIES)
else (WIN32)
	find_library(OPENCL_LIBRARIES OpenCL ENV LD_LIBRARY_PATH)
endif (WIN32)
# if found, add
if (OPENCL_LIBRARIES)
	set(EXTRA_LIBS ${OPENCL_LIBRARIES})
endif (OPENCL_LIBRARIES)

# pointmatcher
find_path(POINT_MATCHER_INCLUDE_DIR pointmatcher/PointMatcher.h
	/usr/local/include
	/usr/include
)
find_library(POINT_MATCHER_LIBRARY libpointmatcher
	/usr/local/lib
	/usr/lib
)

# CGAL
find_path(CGAL_INCLUDE_DIR CGAL/
	/usr/local/include
	/usr/include
)
find_library(CGAL_LIBRARY CGAL
	/usr/local/lib
    /usr/lib
)
if (CGAL_INCLUDE_DIR AND CGAL_LIBRARY)
	set(EXTRA_LIBS ${EXTRA_LIBS} ${CGAL_LIBRARY})
	include_directories(${CGAL_INCLUDE_DIR})
	add_definitions(-DHAVE_CGAL)	
endif (CGAL_INCLUDE_DIR AND CGAL_LIBRARY)
	

# include all libs so far
include_directories(${CMAKE_SOURCE_DIR} ${POINT_MATCHER_INCLUDE_DIR} ${EIGEN2_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${NABO_INCLUDE_DIR})
add_definitions(-Wall)

# main pointmatcher lib
add_library(pointmesher
	pointmesher/Meshers.cpp
	pointmesher/MeshFilters.cpp
)
install(TARGETS pointmesher ARCHIVE DESTINATION lib)
install(FILES
	pointmesher/PointMesher.h
	DESTINATION include/pointmesher
)

# test programs
add_executable(mesh tests/mesh.cpp)
target_link_libraries(mesh pointmesher ${POINT_MATCHER_LIBRARY} ${Boost_LIBRARIES} ${NABO_LIBRARY} ${EXTRA_LIBS})

