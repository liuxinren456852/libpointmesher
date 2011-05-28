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


#include <math.h>

// libpointmatcher & libpointmesher
#include "PointMesher.h"

// Eigenvalues
#include <Eigen/Eigen>

// CGAL
#ifdef HAVE_CGAL
	#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
	//#include <CGAL/Delaunay_triangulation_2.h>
	//#include <CGAL/Triangulation_vertex_base_with_info_2.h>
	//#include <CGAL/centroid.h>
	//#include <CGAL/Simple_cartesian.h>
	//#include <CGAL/Polyhedron_incremental_builder_3.h>
	//#include <CGAL/Polyhedron_3.h>
	//#include <CGAL/Polyhedron_items_with_id_3.h>
	//#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
	//#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
	//#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#endif // HAVE_CGAL


/**
 * MeshFilter base class
 **/

/*
// Compute normal of plane in 3D
template<typename T>
typename PointMesher<T>::Vector3 PointMesher<T>::MeshFilter::computeNormal(Matrix3 matrixIn) const
{
	Vector3 v1 = matrixIn.col(1) - matrixIn.col(0);
	Vector3 v2 = matrixIn.col(2) - matrixIn.col(0);
	Vector3 vn = v1.cross(v2);

	return vn.normalized();
}
*/

/*
#ifdef HAVE_CGAL
// Compute triangle centroid
template<typename T>
typename PointMesher<T>::Vector3 PointMesher<T>::MeshFilter::computeCentroid(const Matrix3 matrixIn) const
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
#endif // HAVE_CGAL
*/


// Compute Euclidean distances of vertices w.r.t. a specified point (default point: origin)
template<typename T>
void PointMesher<T>::MeshFilter::computeVertexDist(const Mesh& meshIn, Vector& vDist, Matrix& mDist, const Vector3 vecPt) const
{
	int nbVert = meshIn.n_vertices();
	vDist.resize(nbVert);
	int nbFaces = meshIn.n_faces();
	mDist.resize(nbFaces, 3);

	MPoint vPoint;
	MPoint vPointOrig(vecPt(0), vecPt(1), vecPt(2));

	// Iterate over all vertices
	for (CVIter vIter = meshIn.vertices_begin(); vIter != meshIn.vertices_end(); ++vIter)
	{
		vPoint = meshIn.point(vIter.handle());
		vDist(vIter.handle().idx()) = (vPoint - vPointOrig).norm();
	}

	// Iterate over all faces
	int indexV = 0;
	for (CFIter fIter = meshIn.faces_begin(); fIter != meshIn.faces_end(); ++fIter)
	{
		// Circulate around the current face
  		for (CFVIter fvIter = meshIn.cfv_iter(fIter.handle()); fvIter; ++fvIter)
		{
			mDist(fIter.handle().idx(), indexV) = vDist(fvIter.handle().idx());
			indexV++;
		}
		indexV = 0;
  	}
}


// Compute perimeters of all the faces
template<typename T>
void PointMesher<T>::MeshFilter::computePerimeters(const Mesh& meshIn, Vector& vPerim, Matrix& mEdges) const
{
	int nbFaces = meshIn.n_faces();
	vPerim.resize(nbFaces);
	mEdges.resize(nbFaces, 3);

	// Iterate over all faces
	VHandle vhTo, vhFrom;
	int indexE = 0;
	for (CFIter fIter = meshIn.faces_begin(); fIter != meshIn.faces_end(); ++fIter)
	{
		// Circulate around the current face
  		for (CFHIter fhIter = meshIn.cfh_iter(fIter.handle()); fhIter; ++fhIter)
		{
			vhTo = meshIn.to_vertex_handle(fhIter.handle());
			vhFrom = meshIn.from_vertex_handle(fhIter.handle());
			mEdges(fIter.handle().idx(), indexE) = (meshIn.point(vhTo) - meshIn.point(vhFrom)).norm();
			indexE++;
		}
		indexE = 0;
	}

	vPerim = mEdges.rowwise().sum();
}


/*
// Compute incident angles of faces w.r.t. a specified point (default: origin)
Vector computeIncAngles(const Vector vPt = Vector3(0, 0, 0)) const
{
	Matrix mCentroids = getFaceAttrByName("Centroids");
	Matrix mNormals = getFaceAttrByName("Normals");

	assert(mCentroids.rows() == vPt.rows());
	Matrix mDir = vPt - mCentroids.colwise();

	Vector vDirNorm = mDir.colwise().norm();
	for (int i = 0; i < ; i++)
	{
		mDir.row(i) = mDir.row(i).cwise() / vDirNorm;
	}

	Vector vIncAngle = (mNormals.cwise() * mDir).colwise().sum;
	vIncAngle = mComp.colwise().sum();
	vIncAngle = vIncAngle.colwise().abs();
	return vIncAngle;
}
*/

template struct PointMesher<float>::MeshFilter;
template struct PointMesher<double>::MeshFilter;



/**
 * Mesh operations
 **/

// Identity
template<typename T>
typename PointMesher<T>::Mesh PointMesher<T>::IdentityMeshFilter::filter(const Mesh& meshIn, bool& iterate)
{
	return meshIn;
}

template struct PointMesher<float>::IdentityMeshFilter;
template struct PointMesher<double>::IdentityMeshFilter;


// Removal of artifacts from mesh

// Constructor
template<typename T>
PointMesher<T>::ArtifactsRemovalMeshFilter::ArtifactsRemovalMeshFilter(
	const T thresh1, const T thresh2, const T thresh3) :
	thresh1(thresh1), thresh2(thresh2), thresh3(thresh3)
{
}
// Filter
template<typename T>
typename PointMesher<T>::Mesh PointMesher<T>::ArtifactsRemovalMeshFilter::filter(const Mesh& meshIn, bool& iterate) const
{
	// Initialization
	Mesh meshOut = meshIn;
	int nbFaces = meshIn.n_faces();
	std::cout << "# faces start: " << nbFaces << std::endl;

	/* Filter 1
	*	for every face in the mesh, apply threshold on the ratio
	*	between its closest and farthest vertex relative to the sensor origin
	*	-> remove shadow faces
	*/
	if (thresh1 > 0)
	{  	
		// Compute distances
		Vector vDist;
		Matrix mDist;
		computeVertexDist(meshOut, vDist, mDist);
		
		Vector mRatio(nbFaces);
		mRatio = mDist.rowwise().maxCoeff();
		mRatio = mRatio.array() / mDist.rowwise().minCoeff().array();
		nbFaces = (mRatio.array() < thresh1).count();

		// Filtering
		for (FIter fIter = meshOut.faces_begin(); fIter != meshOut.faces_end(); ++fIter)
		{
			if (mRatio(fIter.handle().idx()) >= thresh1)
			{
				// Delete face
				meshOut.delete_face(fIter, true);
			}
  		}
	}

	meshOut.garbage_collection();
	nbFaces = meshOut.n_faces();
	std::cout << "# faces filter 1: " << nbFaces << std::endl;
	

	/* Filter 2
	*	for every face in the mesh, apply threshold on the face's perimeter
	*	-> remove frontier faces
	*/

	if (thresh2 > 0 && nbFaces > 0)
	{
		// Compute perimeters
		Vector vPerim;
		Matrix mEdges;
		computePerimeters(meshOut, vPerim, mEdges);
	
		nbFaces = (vPerim.array() < thresh2).count();

		// Filtering
		for (FIter fIter = meshOut.faces_begin(); fIter != meshOut.faces_end(); ++fIter)
		{
			if (vPerim(fIter.handle().idx()) >= thresh2)
			{
				// Delete face
				meshOut.delete_face(fIter, true);
			}
  		}
	}

	meshOut.garbage_collection();
	nbFaces = meshOut.n_faces();
	std::cout << "# faces filter 2: " << nbFaces << std::endl;


	/* Filter 3
	*	for every face in the mesh, apply threshold on the faces
	*	with small incident angle relative to the sensor's line of sight
	* 	-> remove remaining shadow faces that escaped Filter 1
	*/
/*
	if (thresh3 > 0 && j > 0)
	{
		Vector vIncAngle = mesh.computeIncAngles();
		nbF = (vIncAngle.cwise() > thresh3).count();

		// Filtering
		fListOut = Matrix(dimF, nbF);
		fAttrListOut = Matrix(dimAttrF, nbF);

		j = 0;
		for (int i = 0; i < nbF; i++)
		{
			if (vIncAngle(i) > thresh3)
			{
				fListOut.col(j) = fListIn.col(i);
				fAttrListOut.col(j) = fAttrListIn.col(i);
				j++;
			}
		}

		// Update
		fListIn = fListOut;
		fAttrListIn = fAttrListOut;
		Mesh mesh(vListIn, vAttrListIn, meshIn.vertexAttrLabels,
				  fListIn, fAttrListIn, meshIn.faceAttrLabels);
	}
*/	

	std::cout << "# faces end: " << nbFaces << std::endl;

	// Save into vtk-file (for debugging)
	if (!OpenMesh::IO::write_mesh(meshOut, "meshITM_filtered.off")) 
	{
  		std::cerr << "write error\n";
  		exit(1);
	}

	return meshOut;
}

template struct PointMesher<float>::ArtifactsRemovalMeshFilter;
template struct PointMesher<double>::ArtifactsRemovalMeshFilter;


// Simplify mesh
/*
// Modifier creating a surface mesh
typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> SurfMesh;
typedef SurfMesh::HalfedgeDS HalfedgeDS;
typedef CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B;

template <typename TT>
class Build_SM : public CGAL::Modifier_base<HalfedgeDS>
{
	typedef typename Eigen::Matrix<TT, Eigen::Dynamic, Eigen::Dynamic> Matrix;

private:
	int nbV, nbF;
	Matrix mSMPoints;

public:
	Build_SM(Matrix matrixIn)
	{
		mSMPoints = matrixIn;
		nbF = mSMPoints.cols();
		nbV = 3*nbF;//int(mSMPoints.block(12, 0, 3, nbF).maxCoeff() + 1); // upper bound on # vertices
	};
	void operator() (HalfedgeDS & hds)
	{	
		B ib(hds, true);
		ib.begin_surface(nbV, nbF);

		int* arrIndRemap = new int[nbV];
		for (int j = 0; j < nbV; j++)
		{
			arrIndRemap[j] = -1;
		}

		int newIndex, oldIndex, triIndex;
		triIndex = 0;
		for (int i = 0; i < nbF; i++)
		{
			Point p1(mSMPoints(3, i), mSMPoints(4, i), mSMPoints(5, i));
			Point p2(mSMPoints(6, i), mSMPoints(7, i), mSMPoints(8, i));
			Point p3(mSMPoints(9, i), mSMPoints(10, i), mSMPoints(11, i));
			ib.add_vertex(p1);
			ib.add_vertex(p2);
			ib.add_vertex(p3);

			// Rearrange indices
			ib.begin_facet();
			for (int k = 12; k < 15; k++)
			{
				oldIndex = int(mSMPoints(k, i));

				if (arrIndRemap[oldIndex] == -1)
				{
					newIndex = triIndex;
					arrIndRemap[oldIndex] = newIndex;
					triIndex++;
				}
				else
				{
					newIndex = arrIndRemap[oldIndex];
				}

				if (i < 2)
					std::cout << "Remapping: " << oldIndex << " -> " << newIndex << "\n";
				ib.add_vertex_to_facet(newIndex);
				ib.vertex(newIndex)->id() = newIndex;
			}
			ib.end_facet();
		}
		ib.end_surface();
		delete [] arrIndRemap;

		std::cout << triIndex << " " << nbV << "\n";
	}
};
	
// Constructor
template<typename T>
PointMesher<T>::SimplifyMeshingFilter::SimplifyMeshingFilter(
	const int edgeCount) : edgeCount(edgeCount) {}

// Prefilter 
template<typename T>
typename PointMesher<T>::DataPoints PointMesher<T>::SimplifyMeshingFilter::preFilter(
    const DataPoints& input, bool &iterate) const
{
	assert((input.descriptors.rows() >= 15) && (input.features.cols() == input.descriptors.cols()));

	typedef CGAL::Simple_cartesian<double> K;
	typedef K::Point_3 Point;
	typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> SurfMesh;
	typedef SurfMesh::Facet_iterator FacetIt;
	typedef SurfMesh::Halfedge::Halfedge_handle HalfedgeIt;
	typedef SurfMesh::HalfedgeDS HalfedgeDS;
	typedef CGAL::Polyhedron_incremental_builder_3<HalfedgeDS> B;
	namespace SMS = CGAL::Surface_mesh_simplification;

	typename DataPoints::Features featIn;
	typename DataPoints::Descriptors descIn;
	typename DataPoints::Features featOut;
	typename DataPoints::Descriptors descOut;

	featIn = input.features;
	descIn = input.descriptors;

	// Create surface mesh
	SurfMesh sm;
	Build_SM<T> triMesh(descIn);
	sm.delegate(triMesh);

std::cout << sm.size_of_facets() << " Done1\n";

	// Simplify polyhedral mesh
	SMS::Count_stop_predicate<SurfMesh> stop(edgeCount); // stop predicate: stop when number of edges < edgeCount
	int nbEdgesRem;

    std::cout << sm.is_pure_triangle() << " triangle surface mesh.\n";
	std::cout << sm.size_of_facets() << " triangles.\n";


//	nbEdgesRem = SMS::edge_collapse(sm, stop, CGAL::edge_index_map(boost::get(CGAL::edge_external_index, sm)));

std::cout << edgeCount << " " << nbEdgesRem << " Done2\n";

	// Get result
	FacetIt itF; // iterate through facets
	HalfedgeIt itH;
	int j = 0;
	int nbFacets = sm.size_of_facets();
	featOut = typename DataPoints::Features(input.features.rows(), nbFacets);
	descOut = typename DataPoints::Descriptors(input.descriptors.rows(), nbFacets);
	
	Matrix3 vp;
	Vector3 vc;
	Vector3 vn;

	for (itF = sm.facets_begin(); itF != sm.facets_end(); itF++)
	{
		// Access vertices
		itH = itF->halfedge();
		Point p1 = itH->vertex()->point();
		int id1 = itH->vertex()->id();
		Point p2 = itH->next()->vertex()->point();
		int id2 = itH->next()->vertex()->id();
		Point p3 = itH->prev()->vertex()->point();
		int id3 = itH->prev()->vertex()->id();
		vp << p1.x(), p2.x(), p3.x(),
			  p1.y(), p2.y(), p3.y(),
			  p1.z(), p2.z(), p3.z();

		// Compute triangle centroid
		vc = computeCentroid(vp);

		// Compute normal at centroid
		vn = computeNormal(vp);

		featOut.col(j) = vc; 
		descOut.block(0, j, 3, 1) = vn;
		descOut.block(3, j, 3, 1) = vp.col(0);
		descOut.block(6, j, 3, 1) = vp.col(1);
		descOut.block(9, j, 3, 1) = vp.col(2);
		descOut.block(12, j, 3, 1) << id1, id2, id3;
		j++;
	}

std::cout << vp << "\n" << j << " Done3\n";

	return DataPoints(featOut, input.featureLabels, descOut, input.descriptorLabels);
}


template struct PointMesher<float>::SimplifyMeshingFilter;
template struct PointMesher<double>::SimplifyMeshingFilter;
*/

//#endif // HAVE_CGAL
