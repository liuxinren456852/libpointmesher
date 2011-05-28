// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

Generation of triangular surface meshes.


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
#include <vector>
#include <iostream>

// libpointmatcher & libpointmesher
#include "PointMesher.h"

// qhull
#include "qhull/libqhull.h"
//#include "mem.h"
//#include "qset.h"

// PCL
#ifdef HAVE_PCL
	#include <pcl/kdtree/impl/kdtree_flann.hpp>
	#include <pcl/features/normal_3d.h>
	#include <pcl/surface/mls.h>
	#include <pcl/surface/gp3.h>
	#include <pcl/io/vtk_io.h>
	#include <pcl/io/pcd_io.h>
	#include <pcl/filters/statistical_outlier_removal.h>
#endif // HAVE_PCL

// CGAL
#ifdef HAVE_CGAL
	#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
	#include <CGAL/Delaunay_triangulation_2.h>
	#include <CGAL/Triangulation_vertex_base_with_info_2.h>
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
 * Mesher base class
 **/

#ifdef HAVE_PCL
// Conversion from libpointmatcher point cloud to PCL point cloud
template<typename T>
void PointMesher<T>::Mesher::convertPclDatapoints(const DataPoints& ptCloud, DataPointsPCL& ptCloudPCL) const
{
	assert(ptCloud.features.rows() >= 3);
	const int nbPoints = ptCloud.features.cols();
	const int nbDesc = ptCloud.descriptors.cols();

	if (nbPoints == 0)
	{
		std::cout << "[Warning] Mesher: Converted empty point cloud.\n";
	}
	
	std::cout << "# features: " << nbPoints << ", # descriptors: " << nbDesc << std::endl;

	int nbNormals = 0;
	Matrix normals;
	if (nbPoints == nbDesc)
	{
		normals = ptCloud.getDescriptorByName("normals");
		nbNormals = normals.cols();
	}
	
	pcl::PointNormal ptNormal;
	if (nbNormals > 0)
	{
		assert(normals.rows() >= 3);

		// Point cloud consists of points and point normals
		for (int i = 0; i < nbPoints; i++)
		{
			ptNormal.x = ptCloud.features(0, i);
			ptNormal.y = ptCloud.features(1, i);
			ptNormal.z = ptCloud.features(2, i);
			ptNormal.normal[0] = normals(0, i);
			ptNormal.normal[1] = normals(1, i);
			ptNormal.normal[2] = normals(2, i);
			ptCloudPCL.push_back(ptNormal);
		}

		std::cout << "Pts. & normals, # normals: " << nbNormals << std::endl;
	}
	else
	{
		// Point cloud consists of points only
		for (int i = 0; i < nbPoints; i++)
		{
			ptNormal.x = ptCloud.features(0, i);
			ptNormal.y = ptCloud.features(1, i);
			ptNormal.z = ptCloud.features(2, i);
			ptNormal.normal[0] = NAN;
			ptNormal.normal[1] = NAN;
			ptNormal.normal[2] = NAN;
			ptCloudPCL.push_back(ptNormal);
		}
		
		std::cout << "Only Pts., # points: " << nbPoints << std::endl;
	}

	ptCloudPCL.width = nbPoints;
	ptCloudPCL.height = 1;
}


// Conversion from PCL mesh structure to libpointmesher mesh structure
template<typename T>
void PointMesher<T>::Mesher::convertPclPolyMesh(const MeshPCL& triMeshPCL, Mesh& triMesh) const
{
	int nbPoints = triMeshPCL.cloud.width * triMeshPCL.cloud.height;
	int ptSize = triMeshPCL.cloud.data.size() / nbPoints;
    
	int nbFaces = triMeshPCL.polygons.size();
	int nbDim = triMeshPCL.cloud.fields.size();
	float* ptsVal = new float[nbDim];

	typedef typename Mesh::VertexHandle MVHandle;
	MVHandle* vHandle = new MVHandle[nbPoints];
	std::vector<MVHandle> fHandles;

	int vArray[nbPoints];
	for (int k = 0; k < nbPoints; k++)
	{
		vArray[k] = 0;
	}

	int nbPtsPerFace;
	int index;
	for (int i = 0; i < nbFaces; i++)
	{
		nbPtsPerFace = triMeshPCL.polygons[i].vertices.size();
		int* vIndex = new int[nbPtsPerFace];
		for (int j = 0; j < nbPtsPerFace; j++)
		{ 
			index = triMeshPCL.polygons[i].vertices[j];
			vIndex[j] = index;
			if (vArray[index] == 0)
			{
				for (int q = 0; q < nbDim; q++)
				{
					vArray[index] = 1;
					memcpy(ptsVal+q, &triMeshPCL.cloud.data[index*ptSize + triMeshPCL.cloud.fields[q].offset], sizeof(float));
				}
				vHandle[index] = triMesh.add_vertex(MPoint(ptsVal[0], ptsVal[1], ptsVal[2]));

			}
		}

		for (int j = 0; j < nbPtsPerFace; j++)
		{
			fHandles.push_back(vHandle[vIndex[j]]);
		}
		triMesh.add_face(fHandles);

		fHandles.clear();
		delete [] vIndex;
	}

	delete [] ptsVal;
	delete [] vHandle;
}


// Removing outliers using a statistical outlier removal filter
template<typename T>
void PointMesher<T>::Mesher::statOutlierRemovalPCLFilter(const DataPointsPCL& ptCloudPCL_in, DataPointsPCL& ptCloudPCL_out, int mean, double stdMul) const
{
	DataPointsPCL* pPtCloudPCL_in = new DataPointsPCL;
	*pPtCloudPCL_in = ptCloudPCL_in;
	DataPointsPCL::Ptr pPtCloudPCL_shared(pPtCloudPCL_in);
	DataPointsPCL::Ptr pPtCloudPCL_out(new DataPointsPCL);

	// Create the filtering object
	pcl::StatisticalOutlierRemoval<pcl::PointNormal> sor;
	
	// Set filter parameters
	sor.setInputCloud(pPtCloudPCL_shared);
	sor.setMeanK(mean);
	sor.setStddevMulThresh(stdMul);

	// Save negative output (for debugging)
	sor.setNegative(true);
	sor.filter(*pPtCloudPCL_out);
	pcl::PCDWriter writer;
	writer.write<pcl::PointNormal> ("filtertest_outliers.pcd", *pPtCloudPCL_out, false);
	sor.setNegative(false);

	// Filter
	sor.filter(*pPtCloudPCL_out);

	// Save positive output (for debugging)
	writer.write<pcl::PointNormal>("filtertest_inliers.pcd", *pPtCloudPCL_out, false);

	ptCloudPCL_out = *pPtCloudPCL_out;
}


// Smoothing, resampling and normal estimation based on polynomial reconstruction
template<typename T>
void PointMesher<T>::Mesher::mlsResamplingPCLFilter(const DataPointsPCL& ptCloudPCL_in, DataPointsPCL& ptCloudPCL_out, const double searchRadius) const
{
	DataPointsPCL* pPtCloudPCL_in = new DataPointsPCL;
	*pPtCloudPCL_in = ptCloudPCL_in;
	DataPointsPCL::Ptr pPtCloudPCL_shared(pPtCloudPCL_in);
	DataPointsPCL::Ptr pPtCloudPCL_out(new DataPointsPCL);

	// Create KD-tree
	pcl::KdTree<pcl::PointNormal>::Ptr pTree = boost::make_shared<pcl::KdTreeFLANN<pcl::PointNormal> >();
	pTree->setInputCloud(pPtCloudPCL_shared);

	// Create the filtering object
	pcl::MovingLeastSquares<pcl::PointNormal, pcl::PointNormal> mls;
	
	// Set filter parameters
	mls.setOutputNormals(pPtCloudPCL_out);
	mls.setInputCloud(pPtCloudPCL_shared);
	mls.setPolynomialFit(true);
	mls.setSearchMethod(pTree);
	mls.setSearchRadius(searchRadius);

	// Filter
	mls.reconstruct(*pPtCloudPCL_out);

	// Save output (for debugging)
	pcl::io::savePCDFile("filtertest_mls.pcd", *pPtCloudPCL_out);

	ptCloudPCL_out = *pPtCloudPCL_out;
}


// Estimating surface normals in a point cloud using the current viewpoint
template<typename T>
void PointMesher<T>::Mesher::surfaceNormalsPCLFilter(DataPointsPCL& ptCloudPCL, double searchVal, int searchType, double* viewPt) const
{
	// Extract points
	pcl::PointCloud<pcl::PointXYZ>* pPtCloudPCL_points = new pcl::PointCloud<pcl::PointXYZ>;
	pPtCloudPCL_points->points.resize(ptCloudPCL.width);

	for (size_t i = 0; i < pPtCloudPCL_points->size(); ++i)
	{
		pPtCloudPCL_points->points[i].x = ptCloudPCL.points[i].x;
		pPtCloudPCL_points->points[i].y = ptCloudPCL.points[i].y;
		pPtCloudPCL_points->points[i].z = ptCloudPCL.points[i].z;
	}
	pcl::PointCloud<pcl::PointXYZ>::Ptr pPtCloudPCL_shared(pPtCloudPCL_points);

  	// Normal estimation
  	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normEstim;
	pcl::PointCloud<pcl::Normal>::Ptr pPtCloudPCL_normals(new pcl::PointCloud<pcl::Normal>());

  	// Create kdtree representation
  	pcl::KdTree<pcl::PointXYZ>::Ptr pNormTree = boost::make_shared<pcl::KdTreeFLANN<pcl::PointXYZ> >();
	normEstim.setInputCloud(pPtCloudPCL_shared);
	normEstim.setSearchMethod(pNormTree);
	
	// View point
	if (viewPt != NULL)
	{
		normEstim.setViewPoint(viewPt[0], viewPt[1], viewPt[2]);
	}
	else
	{
		normEstim.setViewPoint(0.0, 0.0, 0.0);
	}

	// Type of search 
	if (searchType == 1)
	{
		normEstim.setRadiusSearch(searchVal);
	}
	else // searchType == 2, or default
	{
		normEstim.setKSearch(searchVal);
	}

	normEstim.compute(*pPtCloudPCL_normals);

	// Concatenate the point and normal fields
	pcl::concatenateFields(*pPtCloudPCL_shared, *pPtCloudPCL_normals, ptCloudPCL);

	// Save output (for debugging)
	pcl::io::savePCDFile("filtertest_surfnorm.pcd", ptCloudPCL);
}

template struct PointMesher<float>::Mesher;
template struct PointMesher<double>::Mesher;


/*
template<typename T>
void PointMesher<T>::Mesher::orientNormalsPCLFilter(DataPointsPCL& ptCloud)
{
	// Normal estimation class: flipNormalTowardsViewpoint(...),
	// given the camera acquisition view point (0, 0, 0), do soething like
	// setViewPoint(0, 0, 0) // not necessary, as set to (0,0,0) by default
	// and then estimate the normals towards a specific view point, using:
	// flipNormalTowardsViewpoint(point, 0, 0, 0, normal);
	//
	// Compare to my own code based on libpointmatcher,
	// and more enhanced, come up with spanning tree method
	//
}
*/
#endif // HAVE_PCL


/**
 * Local meshing methods (sensor-centric)
 **/

// Local irregular triangular mesh construction

// Constructor
template<typename T>
PointMesher<T>::ITMLocalMesher::ITMLocalMesher() {}


// Conversions
template<typename T>
typename PointMesher<T>::Matrix PointMesher<T>::ITMLocalMesher::convertCart2Spheric(const Matrix matrixIn) const
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
typename PointMesher<T>::Matrix PointMesher<T>::ITMLocalMesher::computeDelaunay2D(const Matrix matrixIn) const
{
	// ...
	Matrix matrixOut(3, 10);
	
	return matrixOut;
}


// 2D Delaunay triangulation (CGAL)
#ifdef HAVE_CGAL
template<typename T>
typename PointMesher<T>::Matrix PointMesher<T>::ITMLocalMesher::computeDelaunay2D_CGAL(const Matrix matrixIn) const
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
#endif // HAVE_CGAL


// Mesh generation
template<typename T>
typename PointMesher<T>::Mesh PointMesher<T>::ITMLocalMesher::generateMesh(const DataPoints& ptCloud) const
{
	assert(ptCloud.features.rows() >= 3);

	// Get input point cloud
	int nbPoints = ptCloud.features.cols();
	Matrix mVertices = ptCloud.features.block(0, 0, 3, nbPoints);

	// Convert into spherical coordinates
	Matrix mSpheriCoord = convertCart2Spheric(mVertices);

	// 2D Delaunay triangulation
	//Matrix mTriangles = computeDelaunay2D(mSpheriCoord.block(1, 0, 2, mSpheriCoord.cols()));
	Matrix mTriangles = computeDelaunay2D_CGAL(mSpheriCoord.block(1, 0, 2, mSpheriCoord.cols()));

	// Complete ITM by generating a mesh
	Mesh triMesh;
	int nbFaces = mTriangles.cols();
	int nbPtsPerFace = 3;

	VHandle* vHandle = new VHandle[nbPoints];
	std::vector<VHandle> fHandles;

	int index;
	int vIndex[nbPtsPerFace];
	for (int i = 0; i < nbFaces; i++)
	{
		for (int j = 0; j < nbPtsPerFace; j++)
		{ 
			index = mTriangles(j, i);
			vIndex[j] = index;
			vHandle[index] = triMesh.add_vertex(MPoint(mVertices(0, index), mVertices(1, index), mVertices(2, index)));
		}
		for (int j = 0; j < nbPtsPerFace; j++)
		{
			fHandles.push_back(vHandle[vIndex[j]]);
		}
		triMesh.add_face(fHandles);

		fHandles.clear();
	}

	delete [] vHandle;
	
	// Save into vtk-file (for debugging)
	if (!OpenMesh::IO::write_mesh(triMesh, "meshITM.off")) 
	{
  		std::cerr << "write error\n";
  		exit(1);
	}

	return triMesh;
}

template struct PointMesher<float>::ITMLocalMesher;
template struct PointMesher<double>::ITMLocalMesher;


/**
 * Global meshing methods
 **/

// Fast triangulation of unordered point clouds

#ifdef HAVE_PCL
// Constructor
template<typename T>
PointMesher<T>::FastGlobalMesher::FastGlobalMesher(
 		const double searchRadius, const double mu,
		const int maxNN, const double maxSurfAngle,
		const double minAngle, const double maxAngle,
		const bool normConsist, const int knnNormEst) :
		searchRadius(searchRadius), mu(mu),
		maxNN(maxNN), maxSurfAngle(maxSurfAngle),
		minAngle(minAngle), maxAngle(maxAngle),
		normConsist(normConsist), knnNormEst(knnNormEst)
{
}


// Fast triangulation of unordered point clouds
template<typename T>
typename PointMesher<T>::Mesh PointMesher<T>::FastGlobalMesher::generateMesh(const DataPoints& ptCloud) const
{
	// Convert point cloud
	DataPointsPCL* pPtCloudPCL = new DataPointsPCL;
	convertPclDatapoints(ptCloud, *pPtCloudPCL);
	boost::shared_ptr<DataPointsPCL> pPtCloudPCL_shared(pPtCloudPCL);

	// Normal estimation
	DataPointsPCL::Ptr pPtCloudPCL_ptNormals(new DataPointsPCL);
	Matrix normals = ptCloud.getDescriptorByName("normals");
	if (normals.cols() > 0)
	{
		pPtCloudPCL_ptNormals = pPtCloudPCL_shared; //TODO: to be checked if it works like that
	}
	else
	{
		pcl::NormalEstimation<pcl::PointNormal, pcl::Normal> normEstim;
		pcl::PointCloud<pcl::Normal>::Ptr pPtCloudPCL_normals(new pcl::PointCloud<pcl::Normal>());
		pcl::KdTree<pcl::PointNormal>::Ptr pNormTree = boost::make_shared<pcl::KdTreeFLANN<pcl::PointNormal> >();
		pNormTree->setInputCloud(pPtCloudPCL_shared);

		normEstim.setInputCloud(pPtCloudPCL_shared);
		normEstim.setSearchMethod(pNormTree);
		normEstim.setKSearch(knnNormEst);
		normEstim.compute(*pPtCloudPCL_normals);

		// Concatenate the point and normal fields
		pcl::concatenateFields(*pPtCloudPCL_shared, *pPtCloudPCL_normals, *pPtCloudPCL_ptNormals);
	}

	// Create search tree
	pcl::KdTree<pcl::PointNormal>::Ptr pRecTree = boost::make_shared<pcl::KdTreeFLANN<pcl::PointNormal> >();
	pRecTree->setInputCloud(pPtCloudPCL_ptNormals);

	// Initialize objects
	pcl::GreedyProjectionTriangulation<pcl::PointNormal> fastTria;
	pcl::PolygonMesh triMeshPCL;

	// Set the maximum distance between connected points (maximum edge length)
	fastTria.setSearchRadius(searchRadius);

	// Set the other parameters
	fastTria.setMu(mu);
	fastTria.setMaximumNearestNeighbors(maxNN);
	fastTria.setMaximumSurfaceAgle(maxSurfAngle);
	fastTria.setMinimumAngle(minAngle);
	fastTria.setMaximumAngle(maxAngle);
	fastTria.setNormalConsistency(normConsist);

	// Get result
	fastTria.setInputCloud(pPtCloudPCL_ptNormals);
	fastTria.setSearchMethod(pRecTree);
	fastTria.reconstruct(triMeshPCL);

	// Additional vertex information
	//std::vector<int> parts = fastTria.getPartIDs();
	//std::vector<int> states = fastTria.getPointStates();

	// Save into vtk-file (for debugging)
	pcl::io::saveVTKFile("meshFast.vtk", triMeshPCL);

	// Convert mesh
	Mesh triMesh;
	convertPclPolyMesh(triMeshPCL, triMesh);

	return triMesh;
}

template struct PointMesher<float>::FastGlobalMesher;
template struct PointMesher<double>::FastGlobalMesher;

#endif // HAVE_PCL

