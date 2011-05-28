// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

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

#include "PointMesher.h"

#include <math.h>

// Eigenvalues
#include <Eigen/Eigen>

#ifdef HAVE_CGAL
// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/centroid.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>


/**************************************************************************
* Mesh generation
**************************************************************************/

// Mesher

// Compute triangle centroid
template<typename T>
typename PointMesher<T>::Vector3 PointMesher<T>::Mesher::computeCentroid(const Matrix3 matrixIn) const
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
typename PointMesher<T>::Vector3 PointMesher<T>::Mesher::computeNormal(Matrix3 matrixIn) const
{
	Vector3 v1 = matrixIn.col(1) - matrixIn.col(0);
	Vector3 v2 = matrixIn.col(2) - matrixIn.col(0);
	Vector3 vn = v1.cross(v2);

	return vn.normalized();
}

// Create face attributes
template<typename T>
void PointMesher<T>::Mesher::createFaceAttributes(const Matrix vList, const Matrix fList, Matrix& fAttrList, Labels& fAttrLabels) const
{
	// Initialization
	Matrix3 vp;
	Vector3 vc;
	Vector3 vn;

	for (int i = 0; i < fList.cols(); i++)
	{
		// Get triangle
		vp.col(0) = vList.col(int(fList(0, i)));
		vp.col(1) = vList.col(int(fList(1, i)));
		vp.col(2) = vList.col(int(fList(2, i)));

		// Compute triangle centroid
		vc = this->computeCentroid(vp);

		// Compute normal at centroid
		vn = this->computeNormal(vp);

		// Generate data structure
		fAttrList.block(0, i, 3, 1) = vc;
		fAttrList.block(3, i, 3, 1) = vn;
	}

	fAttrLabels[0].text = "FaceCentroids";
	fAttrLabels[0].span = 3;
	fAttrLabels[1].text = "FaceNormals";
	fAttrLabels[1].span = 3;
}

// ---------------------------------

/* Local meshing methods (sensor-centric) */

// Local irregular triangular mesh construction

// Constructor
template<typename T>
PointMesher<T>::ITMLocalMesher::ITMLocalMesher() {}

// Conversions
template<typename T>
typename PointMesher<T>::Matrix PointMesher<T>::ITMLocalMesher::cart2Spheric(const Matrix matrixIn) const
{
	assert(matrixIn.rows() == 3);

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
typename PointMesher<T>::Matrix PointMesher<T>::ITMLocalMesher::delaunay2D(const Matrix matrixIn) const
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

// Mesh generation
template<typename T>
typename PointMesher<T>::Mesh PointMesher<T>::ITMLocalMesher::generateMesh(const DataPoints& ptCloud) const
{
	// Get input point cloud
	dimPts = ptCloud.features.rows();
	nbPts = ptCloud.features.cols();
	Matrix mVertices = ptCloud.features.block(0, 0, (dimPts - 1), nbPts);

	// Convert into spherical coordinates
	Matrix mSpheriCoord = convertCart2Spheric(mVertices);

	// 2D Delaunay triangulation
	Matrix mTriangles = computeDelaunay2D(mSpheriCoord.block(1, 0, 2, mSpheriCoord.cols()));

	// Compute attributes and complete ITM
	Matrix mTriAttr(6, mTriangles.cols());
	Labels mTriLabels(2, Label());
	createFaceAttributes(mVertices, mTriangles, mTriAttr, mTriLabels);

	return Mesh(ptCloud.features, ptCloud.descriptors, ptCloud.descriptorLabels, 
			 	mTriangles, mTriAttr, mTriLabels);
}

template struct PointMesher<float>::ITMLocalMesher;
template struct PointMesher<double>::ITMLocalMesher;

// ---------------------------------

/* Global meshing methods */

// Global irregular triangular mesh construction

// Constructor
template<typename T>
PointMesher<T>::ITMGlobalMesher::ITMGlobalMesher() {}

// Conversions
//...

// Mesh generation
template<typename T>
typename PointMesher<T>::Mesh PointMesher<T>::ITMGlobalMesher::generateMesh(const DataPoints& ptCloud) const
{
	std::cout << "ITMGlobalMesher: Currently identical to ITMLocalMesher!!\n";

	// Get input point cloud
	dimPts = ptCloud.features.rows();
	numPts = ptCloud.features.cols();
	Matrix mVertices = ptCloud.features.block(0, 0, (dimPts - 1), numPts);

	// Convert into spherical coordinates
	Matrix mSpheriCoord = convertCart2Spheric(mVertices);

	// 2D Delaunay triangulation
	Matrix mTriangles = computeDelaunay2D(mSpheriCoord.block(1, 0, 2, mSpheriCoord.cols()));

	// Compute attributes and complete ITM
	Matrix mTriAttr(6, mTriangles.cols());
	Labels mTriLabels(2, Label());
	createFaceAttributes(mVertices, mTriangles, mTriAttr, mTriLabels);

	return Mesh(ptCloud.features, ptCloud.descriptors, ptCloud.descriptorLabels, 
			 	mTriangles, mTriAttr, mTriLabels);
}

template struct PointMesher<float>::ITMLocalMesher;
template struct PointMesher<double>::ITMLocalMesher;

// ---------------------------------

// Global mesh construction by Hoppe's method
//...
	
// Global mesh construction by poisson surface reconstruction
//...

#endif // HAVE_CGAL

