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

#ifndef __POINTMESHER_H
#define __POINTMESHER_H

#include "pointmatcher/PointMatcher.h"

template<typename T>
struct PointMesher
{
	typedef MetricSpaceAligner<T> MSA;
	typedef typename MSA::DataPoints DataPoints;
	typedef typename MSA::Vector Vector;
	typedef typename MSA::Vector3 Vector3;
	typedef typename MSA::Matrix Matrix;
	typedef typename MSA::Matrix3 Matrix3;
	
	// ---------------------------------         

	// Data structure
	struct Mesh
	{
		// Vertices
		Matrix vertexList;
		Matrix vertexAttrList;
		int dimVertex;
		int dimVertexAttr;
		int nbVertices;

		// Faces
		MatrixXi faceList;
		Matrix faceAttrList;
		int dimFace;
		int dimFaceAttr;
		int nbFaces;

		// Labels
		struct Label
		{
			std::string text;
			size_t span;
			Label(const std::string& text = "", const size_t span = 0):text(text), span(span) {}
		};
		typedef std::vector<Label> Labels;

		Labels vertexAttrLabels;
		Labels faceAttrLabels;

		// Constructors
		Mesh() {}

		Mesh(const Matrix& vertexList, const MatriXi& faceList):
			vertexList(vertexList),
			faceList(faceList)
		{}

		Mesh(const Matrix& vertexList, const Matrix& vertexAttrList, const Labels& vertexAttrLabels, 
			 const MatrixXi& faceList, const Matrix& faceAttrList, const Labels& faceAttrLabels):
			vertexList(vertexList),
			vertexAttrList(vertexAttrList),
			vertexAttrLabels(vertexAttrLabels),
			faceList(faceList),
			faceAttrList(faceAttrList),
			faceAttrLabels(faceAttrLabels)
		{}

		/* Getters */
		// Get vertex attribute by name
		// Return a matrix containing only the requested attribute
		Matrix getVertexAttrByName(const std::string& name) const
		{
			int row(0);

			for (unsigned int i = 0; i < vertexAttrLabels.size(); i++)
			{
				const int span(vertexAttrLabels[i].span);
				if (vertexAttrLabels[i].text.compare(name) == 0)
				{
					return vertexAttrList.block(row, 0, span, vertexAttrList.cols());
				}                                                                                      
				row += span;
			}

			return Matrix(); 
		}

		// Get face attribute by name
		// Return a matrix containing only the requested attribute
		Matrix getFaceAttrByName(const std::string& name) const
		{
			int row(0);

			for (unsigned int i = 0; i < faceAttrLabels.size(); i++)
			{
				const int span(faceAttrLabels[i].span);
				if (faceAttrLabels[i].text.compare(name) == 0)
				{
					return faceAttrList.block(row, 0, span, faceAttrList.cols());
				}                                                                                      
				row += span;
			}

			return Matrix(); 
		}

		/* Maintenance functions */
		// Run through list of vertices: remove unconnected vertices
		void cleanVertices()
		{
			// Find unconnected vertices
			int ind, j;
			VectorXi vTag = VectorXi::Zero(nbVertices);
			for (int i = 0; i < nbFaces; i++)
			{
				for (int k = 0; k < dimFace; k++)
				{
					ind = faceList(k, i);
					if (vTag(ind) == 0)
					{
						vTag(ind) = 1;
					}
				}
			}

			// Filter vertices
			vertexListNew = Matrix(dimVertex, nbVertices);
			vertexAttrListNew = Matrix(dimVertexAttr, nbVertices);
			int nbVerticesNew = vTag.sum();

			j = 0;
			for (int q = 0; q < nbVertices; q++)
			{
				if (vTag(q) == 1)
				{
					vertexListNew.col(j) = vertexList.col(q);
					vertexAttrListNew.col(j) = vertexAttrList.col(q);
					j++;
				}
			}

			// Update
			vertexList = vertexListNew;
			vertexAttrList = vertexAttrListNew;
		}
		
		// Run through list of faces: remove invalidated faces, relabel vertices
		void cleanFaces()
		{
			//...
			/*
			for (int i = 0; i < nbFaces; i++)
			{
				for (int k = 0; k < dimFace; k++)
				{
					faceList(k, i) = labelMap(faceList(k, i));
				}
			}
			*/

		}

		/* Additional mesh property functions */
		// Compute perimeter of faces
		Vector computePerimeters() const
		{
			int ind1, ind2;
			Matrix mEdges(nbFaces, dimFace);
			for (int i = 0; i < nbFaces; i++)
			{
				for (int k = 0; k < dimFace; k++)
				{
					ind1 = k;
					ind2 = mod(k+1, dimFace);
					mEdges(i, k) = (vertexList.col(faceList(ind1, i)) - vertexList.col(faceList(ind2, i))).norm();
				}
			}

			Vector vPerim(nbFace);
			vPerim = mEdges.rowwise().sum();
			return vPerim;
		}

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

		// Compute Euclidean distances of vertices w.r.t. a specified point (default: origin)
		Vector computeVertexDist(const Vector vPt = Vector3(0, 0, 0)) const
		{
			Vector vDist(nbVertices);
			vDist = (vertexList.colwise() - vPt).norm();
			return vDist;
		}
	}

	// ---------------------------------         

	// Mesher
	struct Mesher
	{
		virtual ~Mesher() {}
		virtual Mesh generateMesh(const DataPoints& ptCloud) const;

		// TODO: put into a helper class
		Vector3 computeCentroid(const Matrix matrixIn) const;
		Vector3 computeNormal(const Matrix matrixIn) const;

		void createFaceAttributes(const Matrix vList, const Matrix fList, 
								  Matrix& fAttrList, Labels& fAttrLabel) const;
	};

#ifdef HAVE_CGAL

	/* Local meshing methods (sensor-centric) */
	
	// Local irregular triangular mesh construction
	class ITMLocalMesher: public Mesher
	{	
	private:
		// TODO: put into a helper class
		Matrix convertCart2Spheric(const Matrix matrixIn) const;
		Matrix computeDelaunay2D(const Matrix matrixIn) const;
	
	public:
		ITMLocalMesher();
		virtual ~ITMLocalMesher() {};
		virtual Mesh generateMesh(const DataPoints& ptCloud) const;
	};

	/* Global meshing methods */
	
	// Global irregular triangular mesh construction
	class ITMGlobalMesher: public Mesher
	{	
	private:
		// ...
		
	public:
		ITMGlobalMesher();
		virtual ~ITMGlobalMesher() {};
		virtual Mesh generateMesh(const DataPoints& ptCloud) const;
	};

	// Global mesh construction by Hoppe's method
	//class HoppeMesher: public Mesher
	//...
	
	// Global mesh construction by poisson surface reconstruction
	// class PoissonSurfaceMesher: public Mesher
	//...

#endif // HAVE_CGAL

	// ---------------------------------         

	// MeshFilter
	struct MeshFilter
	{
		virtual ~MeshFilter() {}
		virtual void init() {}
		virtual Mesh filter(const Mesh& meshIn, bool& iterate) = 0;

		Vector3 computeCentroid(const Matrix3 matrixIn) const;
	   	Vector3 computeNormal(Matrix3 matrixIn) const;
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
	
#endif // HAVE_CGAL

}; // PointMesher

#endif // __POINTMESHER_H

