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

#include "pointmatcher/PointMatcher.h"
#include "pointmesher/PointMesher.h"

#include <cassert>
#include <iostream>
#include <boost/progress.hpp>


using namespace std;


int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cerr << "Error in command line, usage " << argv[0] << " input_pointcloud.csv (output_mesh.vtk)" << endl;
		return 1;
	}

	typedef MetricSpaceAligner<float> MSA;
	typedef MSA::DataPoints DataPoints;

	typedef PointMesher<float> PM;
	typedef PM::FastGlobalMesher FastGlobalMesher;
	typedef PM::Mesh Mesh;

	DataPoints pc_data = loadCSV<float>(argv[1]);

	/* Datapoints filter */
//	MSA::IdentityDataPointsFilter filter;
	MSA::SamplingSurfaceNormalDataPointsFilter filter(10);
	bool iterate;
	DataPoints subspc_data = filter.filter(pc_data, iterate);

	/* Fast global mesher */
	FastGlobalMesher mesher(0.1, 5, 50, M_PI/2, M_PI/18, 2*M_PI/3, false);
	Mesh* mesh = mesher.generateMesh(subspc_data);
	
	if (argc == 3)
	{
		saveVTK<float>(subspc_data, argv[2]);
	}
	else
	{
		saveVTK<float>(subspc_data, "output_mesh.vtk");
	}

	return 0;
}

