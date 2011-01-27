// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

Copyright (c) 2010--2011,
Andreas Breitenmoserand Stephane Magnenat, ASL, ETHZ, Switzerland
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

#include "PointMesher.h"

#ifdef HAVE_CGAL

/**************************************************************************
* Mesh generation
**************************************************************************/

/** Meshing filter */

// Compute triangle centroid
template<typename T>
typename PointMesher<T>::Vector3 PointMesher<T>::MeshingFilter::computeCentroid(const Matrix3 matrixIn) const
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef K::Point_3 Point;
	typedef CGAL::Triangle_3<K> Triangle;

	// Create triangle	
	Triangle tri(Point(matrixIn(0, 0), matrixIn(1, 0), matrixIn(2, 0)),
				   Point(matrixIn(0, 1), matrixIn(1, 1), matrixIn(2, 1)),
				   Point(matrixIn(0, 2), matrixIn(1, 2), matrixIn(2, 2)));
	Point pointC = centroid(tri);

	Vector3 vc;
	vc(0) = pointC.x();
	vc(1) = pointC.y();
	vc(2) = pointC.z();

	return vc;
}

// Compute normal of plane in 3D
template<typename T>
typename PointMesher<T>::Vector3 PointMesher<T>::MeshingFilter::computeNormal(Matrix3 matrixIn) const
{
	Vector3 v1 = matrixIn.col(1) - matrixIn.col(0);
	Vector3 v2 = matrixIn.col(2) - matrixIn.col(0);
	Vector3 vn = v1.cross(v2);

	return vn.normalized();
}


/** ITMLocalMeshingFilter
  Generate local (i.e. sensor-centric) irregular triangular mesh (ITM) */
 
// Constructor
template<typename T>
PointMesher<T>::ITMLocalMeshingFilter::ITMLocalMeshingFilter() {}

// Conversions
template<typename T>
typename PointMesher<T>::Matrix PointMesher<T>::ITMLocalMeshingFilter::cart2Spheric(const Matrix matrixIn) const
{
	// Generate matrix
	Matrix matrixOut(3, matrixIn.cols());

	// Compute r
	matrixOut.row(0) = matrixIn.colwise().norm();

	for (int i = 0; i < matrixIn.cols(); i++)
	{	
		// Compute theta
		matrixOut(1, i) = asin(matrixIn(2, i) / matrixOut(0, i));
						//acos(matrixIn(2, i) / matrixOut(0, i));

		// Compute phi
		matrixOut(2, i) = atan2(matrixIn(1, i), matrixIn(0, i));
						//acos(matrixIn(0, i) / (matrixOut(0, i)*cos(matrixOut(1, i))));
	}

	return matrixOut;
}

// 2D Delaunay triangulation
template<typename T>
typename PointMesher<T>::Matrix PointMesher<T>::ITMLocalMeshingFilter::delaunay2D(const Matrix matrixIn) const
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> Vb;
	typedef CGAL::Triangulation_face_base_2<K> Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
	typedef CGAL::Delaunay_triangulation_2<K, Tds> DelaunayTri;

	typedef DelaunayTri::Finite_faces_iterator Finite_faces_iterator;
	typedef DelaunayTri::Vertex_handle Vertex_handle;
	typedef DelaunayTri::Point Point;

	// Compute triangulation
	DelaunayTri dt;
  	for (int i = 0; i < matrixIn.cols(); i++)
    {
		// Add point
		Point p = Point(matrixIn(0, i), matrixIn(1, i));
		dt.push_back(Point(p));
		// Add index
		Vertex_handle vh = dt.push_back(p);
		vh->info() = i;
	}

	// Iterate through faces
	Matrix matrixOut(3, dt.number_of_faces());
	int itCol = 0;
	DelaunayTri::Finite_faces_iterator it;
	for (it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++)
	{
		for (int k = 0; k < matrixOut.rows(); k++)
		{
			// Add index of triangle vertex to list
			matrixOut(k, itCol) = it->vertex(k)->info();
		}
		itCol++;
	}

	return matrixOut;
}

// Generate ITM
template<typename T>
void PointMesher<T>::ITMLocalMeshingFilter::generateTriMesh(const Matrix matrixFeatures, const Matrix matrixIndices, 
	Matrix & matrixNewFeatures, Matrix & matrixNewDescriptors) const
{
	// Initialization
	Matrix3 vp;
	Vector3 vc;
	Vector3 vn;

	for (int i = 0; i < matrixIndices.cols(); i++)
	{
		// Get triangle
		vp.col(0) = matrixFeatures.col(int(matrixIndices(0, i)));
		vp.col(1) = matrixFeatures.col(int(matrixIndices(1, i)));
		vp.col(2) = matrixFeatures.col(int(matrixIndices(2, i)));

		// Compute triangle centroid
		vc = computeCentroid(vp);

		// Compute normal at centroid
		vn = computeNormal(vp);

		// Generate data structure
		matrixNewFeatures.col(i) = vc;
		matrixNewDescriptors.block(0, i, 3, 1) = vn;
		matrixNewDescriptors.block(3, i, 3, 1) = vp.col(0);
		matrixNewDescriptors.block(6, i, 3, 1) = vp.col(1);
		matrixNewDescriptors.block(9, i, 3, 1) = vp.col(2);
	}
	matrixNewDescriptors.block(12, 0, 3, matrixNewDescriptors.cols()) = matrixIndices;
}

// Prefilter
template<typename T>
typename PointMesher<T>::DataPoints PointMesher<T>::ITMLocalMeshingFilter::preFilter(
	const DataPoints& input, bool& iterate) const
{
	// Get input features
	Matrix mFeatures = input.features.block(0, 0, 3, input.features.cols());

	// Convert into spherical coordinates
	Matrix mSpheriCoord = cart2Spheric(mFeatures);

	// 2D Delaunay triangulation
	Matrix mTriIndexList = delaunay2D(mSpheriCoord.block(1, 0, 2, mSpheriCoord.cols()));

	// Generate ITM
	Matrix mNewFeatures(3, mTriIndexList.cols());
	Matrix mNewDescriptors(15, mTriIndexList.cols());
	generateTriMesh(mFeatures, mTriIndexList, mNewFeatures, mNewDescriptors);

	return DataPoints(mNewFeatures, input.featureLabels, mNewDescriptors, input.descriptorLabels);
}

template struct PointMesher<float>::ITMLocalMeshingFilter;
template struct PointMesher<double>::ITMLocalMeshingFilter;


/** ITMGlobalMeshingFilter
  Generate global irregular triangular mesh (ITM) */


/** MarchingCubeMeshingFilter
  Generate global surface mesh by Marching Cubes */
  
#endif // HAVE_CGAL
