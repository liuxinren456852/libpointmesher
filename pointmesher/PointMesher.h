// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

Surface reconstruction and geometric modeling library for robotics.


Copyright (c) 2010--2011,
Andreas Breitenmoser and Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the authors at <andreas dot breitenmoser at mavt dot ethz dot ch>
and <stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef __POINTMESHER_H
#define __POINTMESHER_H

// Eigen
#include <Eigen/Eigen>

// Nabo
#include <nabo/nabo.h>

// ICP library --> basic library
#include "pointmatcher/PointMatcher.h"

// OpenMesh & IsoEx
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/VectorT.hh>
//#include <Implicits/ImplicitSphere.hh>
//#include <Implicits/CSG.hh>
//#include <Grids/ImplicitGrid.hh>

// PCL
#ifdef HAVE_PCL
	#include <pcl/point_cloud.h>
	#include <pcl/point_types.h>
	#include <pcl/PolygonMesh.h>
#endif // HAVE_PCL

// CGAL
#ifdef HAVE_CGAL
	//#include ...
#endif // HAVE_CGAL


template<typename T>
struct PointMesher
{

	// Define traits
	struct MeshTraits : public OpenMesh::DefaultTraits
	{
		typedef OpenMesh::VectorT<T, 3> Point;
		typedef OpenMesh::VectorT<T, 3> Normal;
		typedef OpenMesh::VectorT<T, 2> TexCoord;

		VertexAttributes(OpenMesh::Attributes::Status);
  		FaceAttributes(OpenMesh::Attributes::Status);
  		EdgeAttributes(OpenMesh::Attributes::Status);
	};

	typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> Mesh;
	typedef typename Mesh::Point MPoint;
	typedef typename Mesh::VertexHandle VHandle;	
	
	// Iterators & circulators
	typedef typename Mesh::VertexIter VIter;
	typedef typename Mesh::ConstVertexIter CVIter;
	typedef typename Mesh::FaceIter FIter;
	typedef typename Mesh::ConstFaceIter CFIter;
	typedef typename Mesh::FaceVertexIter FVIter;
	typedef typename Mesh::ConstFaceVertexIter CFVIter;
	typedef typename Mesh::FaceHalfedgeIter FHIter;
	typedef typename Mesh::ConstFaceHalfedgeIter CFHIter;


	#ifdef HAVE_PCL
		typedef pcl::PointCloud<pcl::PointNormal> DataPointsPCL;
		typedef pcl::PolygonMesh MeshPCL;
	#endif // HAVE_PCL

	typedef MetricSpaceAligner<T> MSA;
    typedef typename MSA::DataPoints DataPoints;
	typedef typename MSA::Matrix Matrix;

	//typedef Mesh::Scalar MeshScalar;
	//typedef OpenMesh::Vec3d MeshVec3d;
	//typedef typename OpenMesh::VectorT<T, 3> vPoint;

	typedef typename MSA::Vector Vector;
	typedef typename MSA::Vector3 Vector3;
	typedef typename MSA::Matrix3 Matrix3;
	

/**********************************************************************************
 * Mesher
 **********************************************************************************/

	struct Mesher
	{
		virtual ~Mesher() {}
		virtual Mesh generateMesh(const DataPoints& ptCloud) const = 0; 

		#ifdef HAVE_PCL
			/* Conversion of data structures between PCL and libpointmatcher/mesher */
			void convertPclDatapoints(const DataPoints& ptCloud, DataPointsPCL& ptCloudPCL) const; // conversion from DataPoints to PCL PointCloud
			void convertPclPolyMesh(const MeshPCL& triMeshPCL, Mesh& triMesh) const; // conversion from PCL PolygonMesh to Mesh
			
			/* Additional standard point cloud filters from PCL */
			void statOutlierRemovalPCLFilter(const DataPointsPCL& ptCloudPCL_in, DataPointsPCL& ptCloudPCL_out, int mean, double stdMul) const; // wrapper for statistical outlier removal
			void mlsResamplingPCLFilter(const DataPointsPCL& ptCloudPCL_in, DataPointsPCL& ptCloudPCL_out, const double searchRadius) const; // wrapper for Robust Moving Least Squares (RMLS)
			void surfaceNormalsPCLFilter(DataPointsPCL& ptCloud, double searchVal = 10, int searchType = 2, double* viewPt = 0) const; // wrapper for surface normal estimation
			//void orientNormalsPCLFilter(DataPointsPCL& ptCLoud); // wrapper for orientation of surface normals by view point
		#endif // HAVE_PCL
	};


	/**
	 * Local meshing methods (sensor-centric)
	 **/

	// Local irregular triangular mesh construction
	class ITMLocalMesher: public Mesher
	{
		private:
			Matrix convertCart2Spheric(const Matrix matrixIn) const;
			Matrix computeDelaunay2D(const Matrix matrixIn) const;
			#ifdef HAVE_CGAL
				Matrix computeDelaunay2D_CGAL(const Matrix matrixIn) const;
			#endif // HAVE_CGAL

		public:
			ITMLocalMesher();
			virtual ~ITMLocalMesher() {};
			virtual Mesh generateMesh(const DataPoints& ptCloud) const;
	};


	/**
	 * Global meshing methods
	 **/

	#ifdef HAVE_PCL
	// Fast triangulation of unordered point clouds

	class FastGlobalMesher: public Mesher
	{
	private:
		const double searchRadius;
		const double mu;
		const int maxNN;
		const double maxSurfAngle;
		const double minAngle;
		const double maxAngle;
		const bool normConsist;
		const int knnNormEst;

	public:
		// maxSurfAngle = 45 deg, minAngle = 10 deg, maxAngle = 120 deg
		FastGlobalMesher(
			const double searchRadius = 0.025,
			const double mu = 2.5,
			const int maxNN = 100,
			const double maxSurfAngle = M_PI/4,
			const double minAngle = M_PI/18,
			const double maxAngle = 2*M_PI/3,
			const bool normConsist = false,
			const int knnNormEst = 20);
		virtual ~FastGlobalMesher() {};
		virtual Mesh generateMesh(const DataPoints& ptCloud) const;
	};
	#endif // HAVE_PCL


/**********************************************************************************
 * MeshFilter
 **********************************************************************************/

	struct MeshFilter
	{
		virtual ~MeshFilter() {}
		virtual void init() {}
		virtual Mesh filter(const Mesh& meshIn, bool& iterate) const = 0;
		
		/* Computation of mesh properties */
		//Vector3 computeNormal(Matrix3 matrixIn) const; // 
		//#ifdef HAVE_CGAL
			//Vector3 computeCentroid(const Matrix3 matrixIn) const;
		//#endif // HAVE_CGAL
		void computeVertexDist(const Mesh& meshIn, Vector& vDist, Matrix& mDist, const Vector3 vecPt = Vector3(0.0, 0.0, 0.0)) const; // Euclidean distances of vertices w.r.t. a specified point (default point: origin)
		void computePerimeters(const Mesh& meshIn, Vector& vPerim, Matrix& mEdges) const; // perimeter of faces
		//Vector computeIncAngles(const Vector3 vPt = Vector3(0.0, 0.0, 0.0)) const; // incident angles of faces w.r.t. a specified point (default point: origin)
	
	};
	
	struct MeshFilters: public std::vector<MeshFilter*>
	{
		void init();
		void apply(Mesh& mesh, bool iterate);
	};
	typedef typename MeshFilters::iterator MeshFiltersIt;
	typedef typename MeshFilters::const_iterator MeshFiltersConstIt;


	/* Mesh operations */
	
	// Identity
	struct IdentityMeshFilter: public MeshFilter
	{
		virtual Mesh filter(const Mesh& meshIn, bool& iterate);
	};


	#ifdef HAVE_CGAL
	// Removal of artifacts from mesh
	class ArtifactsRemovalMeshFilter: public MeshFilter 
	{
		// Filter thresholds
		const T thresh1;
		const T thresh2;
		const T thresh3;

	public:
		ArtifactsRemovalMeshFilter(const T thresh1 = 1.1, const T thresh2 = 10, const T thresh3 = 0.2);
		virtual ~ArtifactsRemovalMeshFilter() {};
		virtual Mesh filter(const Mesh& meshIn, bool& iterate) const;
 	};

/*
	// Simplify mesh
	class SimplifyMeshFilter: public MeshFilter
	{
		// Filter thresholds
		const int edgeCount;

	public:
		SimplifyMeshFilter(const int edgeCount = 1000);
		virtual ~SimplifyMeshFilter() {};
		virtual Mesh filter(const Mesh& meshIn, bool& iterate) const;
	};
*/
	#endif // HAVE_CGAL

}; // PointMesher

#endif // __POINTMESHER_H

