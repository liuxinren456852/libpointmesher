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


// Compute triangle centroid
template<typename T>
typename PointMesher<T>::Vector3 PointMesher<T>::computeCentroid(const Matrix3 matrixIn) const
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

#endif // HAVE_CGAL

