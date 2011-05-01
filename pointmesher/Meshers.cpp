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

#include "PointMesher.h"

// Eigen
#include <Eigen/Eigen>

// PCL
#ifdef HAVE_PCL
	#include <pcl/point_types.h>
	#include <pcl/kdtree/impl/kdtree_flann.hpp>
	#include <pcl/features/normal_3d.h>
	#include <pcl/surface/mls.h>
	#include <pcl/surface/gp3.h>
	#include <pcl/io/vtk_io.h>
	#include <pcl/io/pcd_io.h>
	#include <pcl/filters/statistical_outlier_removal.h>
#endif

// CGAL
#ifdef HAVE_CGAL
	//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
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
#endif


/**
 * Mesher base class
 **/

#ifdef HAVE_PCL
// Conversion from libpointmatcher point cloud to PCL point cloud.
template<typename T>
typename PointMesher<T>::DataPointsPCL* PointMesher<T>::Mesher::convertPclDatapoints(const DataPoints& ptCloud)
{
	DataPointsPCL* ptCloudConv_ptr = new DataPointsPCL;

	const int nbPoints = ptCloud.features.cols();
	const int nbDesc = ptCloud.descriptors.cols();

	if (nbPoints == 0)
	{
		std::cout << "[Warning] Mesher: Converted empty point cloud.\n";
		return ptCloudConv_ptr;
	}
	
	//std::cout << "# points: " << nbDesc << std::endl;

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
		// Point cloud consists of points and point normals
		int i;
		for (i = 0; i < nbPoints; i++)
		{
			ptNormal.x = ptCloud.features(0, i);
			ptNormal.y = ptCloud.features(1, i);
			ptNormal.z = ptCloud.features(2, i);
			ptNormal.normal[0] = normals(0, i);
			ptNormal.normal[1] = normals(1, i);
			ptNormal.normal[2] = normals(2, i);
			ptCloudConv_ptr->push_back(ptNormal);
		}

		//std::cout << "Pts. & normals, # normals: " << nbNormals << " " << i << std::endl;
	}
	else
	{
		int i;
		// Point cloud consists of points only
		for (i = 0; i < nbPoints; i++)
		{
			ptNormal.x = ptCloud.features(0, i);
			ptNormal.y = ptCloud.features(1, i);
			ptNormal.z = ptCloud.features(2, i);
			//ptNormal.normal[0] = 0.0;
			//ptNormal.normal[1] = 0.0;
			//ptNormal.normal[2] = 0.0;
			ptCloudConv_ptr->push_back(ptNormal);
		}
		
		//std::cout << "Only Pts., # points: " << nbNormals << " " <<  i << std::endl;
	}
	
	return ptCloudConv_ptr;
}


// Conversion from PCL mesh structure to libpointmesher mesh structure.
template<typename T>
typename PointMesher<T>::Mesh* PointMesher<T>::Mesher::convertPclPolyMesh(const MeshPCL& triMeshPCL)
{
	int nbPoints = triMeshPCL.cloud.width * triMeshPCL.cloud.height;
	int ptSize = triMeshPCL.cloud.data.size() / nbPoints;
    
	int nbFaces = triMeshPCL.polygons.size();
	int nbDim = triMeshPCL.cloud.fields.size();
	float* ptsVal = new float[nbDim];

	Mesh triMesh;
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

	return &triMesh;
}


// Removing outliers using a statistical outlier removal filter
template<typename T>
void PointMesher<T>::Mesher::statOutlierRemovalPCLFilter(DataPointsPCL& ptCloudPCL, int mean, double stdMul)
{
	boost::shared_ptr<DataPointsPCL> ptCloudPCL_shared(&ptCloudPCL);
	DataPointsPCL::Ptr ptCloudPCL_filtered(new DataPointsPCL);

	// Create the filtering object
	pcl::StatisticalOutlierRemoval<pcl::PointNormal> sor;
	sor.setInputCloud(ptCloudPCL_shared);
	sor.setMeanK(mean);
	sor.setStddevMulThresh(stdMul);
	sor.filter(*ptCloudPCL_filtered);

	pcl::PCDWriter writer;
	writer.write<pcl::PointNormal>("table_scene_lms400_inliers.pcd", *ptCloudPCL_filtered, false);

	sor.setNegative(true);
	sor.filter(*ptCloudPCL_filtered);
	writer.write<pcl::PointNormal> ("table_scene_lms400_outliers.pcd", *ptCloudPCL_filtered, false);
}


// Smoothing and normal estimation based on polynomial reconstruction
template<typename T>
void PointMesher<T>::Mesher::mlsResamplingPCLFilter(DataPointsPCL& ptCloudPCL, double searchRadius)
{
	boost::shared_ptr<DataPointsPCL> ptCloudPCL_shared(&ptCloudPCL);
	DataPointsPCL::Ptr ptCloudPCL_filtered(new DataPointsPCL);

	// Create KD-tree
	pcl::KdTree<pcl::PointNormal>::Ptr tree = boost::make_shared<pcl::KdTreeFLANN<pcl::PointNormal> >();
	tree->setInputCloud(ptCloudPCL_shared);

	// Create the filtering object
	pcl::MovingLeastSquares<pcl::PointNormal, pcl::PointNormal> mls;
	
	// Set filter parameters
	mls.setOutputNormals(ptCloudPCL_filtered);
	mls.setInputCloud(ptCloudPCL_shared);
	mls.setPolynomialFit(true);
	mls.setSearchMethod(tree);
	mls.setSearchRadius(searchRadius);

	// Filter
	mls.reconstruct(*ptCloudPCL_filtered);

	// Save output
	pcl::io::savePCDFile("bun0-mls.pcd", *ptCloudPCL_filtered);
}

template<typename T>
void PointMesher<T>::Mesher::surfaceNormalsPCLFilter(DataPointsPCL& ptCloud)
{
	// above or first part of fast surface reconstruction
	// Normal estimation
	/*
	NormalEstimation<PointXYZ, Normal> n;
	PointCloud<Normal>::Ptr normals (new PointCloud<Normal>);
	KdTree<PointXYZ>::Ptr tree (new KdTreeFLANN<PointXYZ>);
	tree->setInputCloud (cloud);
	n.setInputCloud (cloud);
	n.setSearchMethod (tree);
	n.setKSearch (20);
	n.compute (*normals);

	// Concatenate the XYZ and normal fields
	PointCloud<PointNormal>::Ptr cloud_with_normals (new PointCloud<PointNormal>);
	pcl::concatenateFields (*cloud, *normals, *cloud_with_normals);
*/
}

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

#endif // HAVE_PCL


/* Local meshing methods (sensor-centric) */

// Local irregular triangular mesh construction
/*
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
*/


/**
 * Global meshing methods
 **/

// Fast triangulation of unordered point clouds

#ifdef HAVE_PCL
// Constructor
template<typename T>
PointMesher<T>::FastGlobalMesher::FastGlobalMesher(const double searchRadius,
		const double mu, const int maxNN,
		const double maxSurfAngle, const double minAngle,
		const double maxAngle, const bool normConsist) :
		searchRadius(searchRadius), mu(mu), maxNN(maxNN),
		maxSurfAngle(maxSurfAngle), minAngle(minAngle),
		maxAngle(maxAngle), normConsist(normConsist)
{
}

// Fast triangulation of unordered point clouds
template<typename T>
typename PointMesher<T>::Mesh* PointMesher<T>::FastGlobalMesher::generateMesh(const DataPoints& ptCloud)
{
	pcl::PointCloud<pcl::PointNormal>* cloud_pcl = new pcl::PointCloud<pcl::PointNormal>();
	//pcl::PointCloud<pcl::PointNormal>::Ptr cloud_pcl(new pcl::PointCloud<pcl::PointNormal>);
	cloud_pcl = Mesher::convertPclDatapoints(ptCloud);
	boost::shared_ptr<pcl::PointCloud<pcl::PointNormal> > cloud_sharedPCL(cloud_pcl);

	// Normal estimation
	pcl::NormalEstimation<pcl::PointNormal, pcl::Normal> normEstim;
	pcl::PointCloud<pcl::Normal>::Ptr cloud_normals(new pcl::PointCloud<pcl::Normal>);
	pcl::KdTree<pcl::PointNormal>::Ptr normTree = boost::make_shared<pcl::KdTreeFLANN<pcl::PointNormal> >();
	normTree->setInputCloud(cloud_sharedPCL);

	normEstim.setInputCloud(cloud_sharedPCL);
	normEstim.setSearchMethod(normTree);
	normEstim.setKSearch(20);
	normEstim.compute(*cloud_normals);

	// Concatenate the point and normal fields
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_ptsNormals(new pcl::PointCloud<pcl::PointNormal>);
	pcl::concatenateFields(*cloud_sharedPCL, *cloud_normals, *cloud_ptsNormals);

	// Create search tree
	pcl::KdTree<pcl::PointNormal>::Ptr recTree = boost::make_shared<pcl::KdTreeFLANN<pcl::PointNormal> >();
	recTree->setInputCloud(cloud_ptsNormals);

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
	fastTria.setInputCloud(cloud_ptsNormals);
	fastTria.setSearchMethod(recTree);
	fastTria.reconstruct(triMeshPCL);

	// Additional vertex information
	std::vector<int> parts = fastTria.getPartIDs();
	std::vector<int> states = fastTria.getPointStates();

std::cout << "Hello1" << std::endl;

	pcl::io::saveVTKFile("mesh.vtk", triMeshPCL);

std::cout << "Hello2" << std::endl;

	Mesh* triMesh_ptr;

std::cout << "Hello3" << std::endl;

	triMesh_ptr = Mesher::convertPclPolyMesh(triMeshPCL);

std::cout << "Hello4" << std::endl;

	return triMesh_ptr; 
}

template struct PointMesher<float>::FastGlobalMesher;
template struct PointMesher<double>::FastGlobalMesher;

#endif // HAVE_PCL

