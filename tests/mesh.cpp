// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab

#include "pointmatcher/PointMatcher.h"
#include "pointmesher/PointMesher.h"
#include <cassert>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{

	#ifdef HAVE_CGAL

	if (argc != 2)
	{
		cerr << "Error in command line, usage " << argv[0] << " cloud.csv" << endl;
		return 1;
	}
		
	typedef MetricSpaceAlignerD MSA;
	typedef MSA::DataPoints DP;
	
	MSA::VTKFileInspector inspector("output");

	/* Load data */
	DP data(loadCSV<MSA::ScalarType>(argv[1]));
	inspector.dumpDataPoints(data, "point_cloud"); // dump data 1

	/* Point cloud filters */
	MSA::DataPointsFilters filtersCloud;
	bool iterate(true);
	//filtersCloud.push_back(new MSA::FixstepSamplingDataPointsFilter);
	filtersCloud.push_back(new MSA::FixstepSamplingDataPointsFilter(100));
	filtersCloud.push_back(new MSA::SurfaceNormalDataPointsFilter(10, 0));
	filtersCloud.applyPre(data, iterate);
	inspector.dumpDataPoints(data, "point_cloud_subsampled"); // dump data 2

	/* Meshing filters */
	MSA::DataPointsFilters filtersMesh;
	filtersMesh.push_back(new MSA::ITMLocalMeshingFilter);
	filtersMesh.push_back(new MSA::ArtifactsRemovalMeshingFilter(1.1, 0.5, 0.2));
	//filtersMesh.push_back(new MSA::ArtifactsRemovalMeshingFilter(1.25, 10, 0.15));
	//filtersMesh.push_back(new MSA::SimplifyMeshingFilter(10000));
	filtersMesh.applyPre(data, iterate);
	inspector.dumpDataPoints(data, "mesh_centroid");
	inspector.dumpMeshNodes(data, "triangle_mesh"); // dump data 3

	#endif // HAVE_CGAL

	return 0;
}

