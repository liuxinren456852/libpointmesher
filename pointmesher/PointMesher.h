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
	
	// TODO: separate into Mesher and MeshFilter
	// if needed, add helper class/functions/interfaces
	
	#ifdef HAVE_CGAL
	
	// MeshingFilter
	struct MeshingFilter: public MSA::DataPointsFilter
	{
		Vector3 computeCentroid(const Matrix3 matrixIn) const;
		Vector3 computeNormal(Matrix3 matrixIn) const;
		
		virtual ~MeshingFilter() {};
		virtual DataPoints preFilter(const DataPoints& input, bool& iterate) const = 0;
		virtual DataPoints stepFilter(const DataPoints& input, bool& iterate) const = 0;
	};

	// ITMLocalMeshingFilter
	class ITMLocalMeshingFilter: public MeshingFilter
	{
	public:
		ITMLocalMeshingFilter();
		virtual ~ITMLocalMeshingFilter() {};
		virtual DataPoints preFilter(const DataPoints& input, bool& iterate) const;
		virtual DataPoints stepFilter(const DataPoints& input, bool& iterate) const {return input;};

	private:
		Matrix cart2Spheric(const Matrix matrixIn) const;
		Matrix delaunay2D(const Matrix matrixIn) const;
		void generateTriMesh(const Matrix matrixFeatures, const Matrix matrixIndices, Matrix& matrixNewFeatures, Matrix& matrixNewDescriptors) const;
	};

	// ITMGlobalMeshingFilter
	// ...

	// MarchingCubeMeshingFilter
	// ...

	// ArtifactsRemovalMeshingFilter
	class ArtifactsRemovalMeshingFilter: public MeshingFilter 
	{
		// Filter thrasholds
		const T thresh1;
		const T thresh2;
		const T thresh3;

	public:
		ArtifactsRemovalMeshingFilter(const T thresh1 = 1.1, const T thresh2 = 10, const T thresh3 = 0.2);
		virtual ~ArtifactsRemovalMeshingFilter() {};
		virtual DataPoints preFilter(const DataPoints& input, bool& iterate) const;
		virtual DataPoints stepFilter(const DataPoints& input, bool& iterate) const {return input;};
 	};

	// SimplifyMeshingFilter
	class SimplifyMeshingFilter: public MeshingFilter
	{
		// Filter thrasholds
		const int edgeCount;

		public:
		SimplifyMeshingFilter(const int edgeCount = 1000);
		virtual ~SimplifyMeshingFilter() {};
		virtual DataPoints preFilter(const DataPoints& input, bool& iterate) const;
		virtual DataPoints stepFilter(const DataPoints& input, bool& iterate) const {return input;};
	};
	
	#endif // HAVE_CGAL
};

#endif // __POINTMESHER_H
