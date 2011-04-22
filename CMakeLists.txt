cmake_minimum_required(VERSION 2.6)

project(libpointmesher)

#include(UseDoxygen)


list(APPEND CMAKE_MODULE_PATH ${libpointmesher_SOURCE_DIR}/cmake)


# Boost library
#########################################
find_package(Boost REQUIRED)


# Linear algebra library
#########################################

# Eigen 3
find_path(EIGEN_INCLUDE_DIR Eigen/Core
			/usr/local/include/eigen
			/usr/include/eigen
)

## Eigen 2
#find_path(EIGEN2_INCLUDE_DIR Eigen/Core
#			/usr/local/include/eigen2
#			/usr/include/eigen2
#)


# Nearest neighbor search library
#########################################

# Nabo
find_path(NABO_INCLUDE_DIR nabo/nabo.h
			/usr/local/include
			/usr/include
)
find_library(NABO_LIBRARY libnabo
			/usr/local/lib
			/usr/lib
)

# Flann
find_path(FLANN_INCLUDE_DIR src/cpp
 			/usr/local/include/flann
 			/usr/include/flann
)
find_library(FLANN_LIBRARY libflann_cpp_s
			/usr/local/lib
 			/usr/lib
)
if (FLANN_INCLUDE_DIR AND FLANN_LIBRARY)
	set(PCL_LIBS ${PCL_LIBS} ${FLANN_LIBRARY})
	include_directories(${FLANN_INCLUDE_DIR})
	set(FLANN_OK true)
endif (FLANN_INCLUDE_DIR AND FLANN_LIBRARY)


# ICP library
#########################################

# Pointmatcher
find_path(POINT_MATCHER_INCLUDE_DIR pointmatcher/PointMatcher.h
			/usr/local/include
			/usr/include
)
find_library(POINT_MATCHER_LIBRARY libpointmatcher
			/usr/local/lib
			/usr/lib
)


# Equation solver library
#########################################

# CMinpack 
find_path(CMINPACK_INCLUDE_DIR .
			/usr/local/include/cminpack   
    		/usr/include/cminpack
)
if (CMINPACK_INCLUDE_DIR)
    set(PCL_LIBS ${PCL_LIBS} ${CMINPACK_LIBRARY})
	include_directories(${CMINPACK_INCLUDE_DIR})
	set(CMINPACK_OK true)
endif (CMINPACK_INCLUDE_DIR)


# Computational geometry library
#########################################

# Qhull
find_path(QHULL_INCLUDE_DIR src2
			/usr/local/include/qhull
			/usr/include/qhull
)
find_library(QHULL_LIBRARY libqhull
			/usr/local/lib
			/usr/lib
)
if (QHULL_INCLUDE_DIR AND QHULL_LIBRARY)
    set(PCL_LIBS ${PCL_LIBS} ${QHULL_LIBRARY})
	include_directories(${QHULL_INCLUDE_DIR})
	set(QHULL_OK true)
endif (QHULL_INCLUDE_DIR AND QHULL_LIBRARY)


# Computer graphics library
#########################################

# OpenMesh & IsoEx
find_path(OPENMESH_INCLUDE_DIR OpenMesh
    		/usr/local/include/openmesh
	    	/usr/include/openmesh
)
find_library(OPENMESH_LIBRARY libOpenMeshCore
		    /usr/local/lib
			/usr/lib
)
find_path(ISOEX_INCLUDE_DIR src/Math
							src/Implicits
							src/Grids
							src/Extractors
							/usr/local/include/isoex
							/usr/include/isoex
)
find_library(ISOEX_LIBRARY libisoex
			/usr/local/lib
			/usr/lib
)

# PCL
if (FLANN_OK AND CMINPACK_OK AND QHULL_OK)
	find_package(PCL 1.0 REQUIRED)
	if (PCL_FOUND)
		set(EXTRA_LIBS ${PCL_LIBS} ${PCL_IO_LIBRARY}
								 ${PCL_FILTERS_LIBRARY}
								 ${PCL_FEATURES_LIBRARY}
								 ${PCL_SEGMENTATION_LIBRARY}
								 ${PCL_SURFACE_LIBRARY}
								 ${PCL_REGISTRATION_LIBRARY}
#								 ${PCL_KEYPOINTS_LIBRARY}
#								 ${PCL_RANGEIMAGE_LIBRARY}
		   )
		include_directories(${PCL_INCLUDE_DIRS})
		set(add_definitions(-DHAVE_PCL))
	endif (PCL_FOUND)
endif (FLANN_OK AND CMINPACK_OK AND QHULL_OK)

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


# OpenCL library
#########################################

if (WIN32)
	find_library(OPENCL_LIBRARIES opencl64)
	if (!OPENCL_LIBRARIES)
		find_library(OPENCL_LIBRARIES opencl32)
	endif (!OPENCL_LIBRARIES)
else (WIN32)
	find_library(OPENCL_LIBRARIES OpenCL ENV LD_LIBRARY_PATH)
endif (WIN32)
if (OPENCL_LIBRARIES)
	set(EXTRA_LIBS ${EXTRA_LIBS} ${OPENCL_LIBRARIES})
endif (OPENCL_LIBRARIES)


# Include all libraries so far
#########################################

include_directories(${CMAKE_SOURCE_DIR} ${EIGEN_INCLUDE_DIR} ${NABO_INCLUDE_DIR} ${POINT_MATCHER_INCLUDE_DIR} ${BOOST_INCLUDE_DIRS} ${OPENMESH_INCLUDE_DIR} ${ISOEX_INCLUDE_DIR})
#${EIGEN2_INCLUDE_DIR})
add_definitions(-Wall)


# Main pointmatcher library
#########################################

add_library(pointmesher
	pointmesher/Meshers.cpp
#	pointmesher/MeshFilters.cpp
#	pointmesher/MeshHelpers.cpp
)
install(TARGETS pointmesher ARCHIVE DESTINATION lib)
install(FILES
	pointmesher/PointMesher.h
	DESTINATION include/pointmesher
)


# Test programs
#########################################

add_executable(mesh tests/mesh.cpp)
target_link_libraries(mesh pointmesher ${POINT_MATCHER_LIBRARY} ${NABO_LIBRARY} ${BOOST_LIBRARIES} ${OPENMESH_LIBRARY} ${ISOEX_LIBRARY} ${EXTRA_LIBS})

