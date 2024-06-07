#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_cloud_io.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"

#include "PoissonReconstruct.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_grid.h"

using namespace geometrycentral;
using namespace geometrycentral::pointcloud;



static void convertToEigenDS(const Vector<Vector3>& in, Eigen::MatrixXd& out)
{
    out.resize(in.rows(), 3);
    for (size_t i = 0; i < in.rows(); ++i)
    {
        out(i, 0) = in[i].x;
        out(i, 1) = in[i].y;
        out(i, 2) = in[i].z;
    }
}

int main() {
    std::unique_ptr<PointCloud> cloud;
    std::unique_ptr<PointPositionGeometry> geom;
    std::string path = "D:/dev/GeometricAlgo/input/PointClouds/kitten_04.obj";
    std::tie(cloud, geom) = readPointCloud(path);

    PointCloud* cloudP = cloud.release();
    PointPositionGeometry* geoP = geom.release();


    polyscope::init();

    // Register a point cloud
    // `points` is a Nx3 array-like container of points
    polyscope::registerPointCloud("my points", geoP->positions);

    geoP->requireNormals();
    // prepare for poisson-reconstruct
    Eigen::MatrixXd P, N;
    Eigen::Vector3i dims;
    std::pair<Eigen::Vector3d, Eigen::Vector3d> bbox;
    double iso;
    Eigen::MatrixXd x;
    Eigen::VectorXd g;

    auto rawPoints = geoP->positions.raw();
    auto rawNormals = -geoP->normals.raw();
    convertToEigenDS(rawPoints, P);
    convertToEigenDS(rawNormals, N);
    PoissonRec::poisson_rec(P, N, dims, bbox, iso, x, g);
    uint32_t dimX = dims.x();
    uint32_t dimY = dims.y();
    uint32_t dimZ = dims.z();
    glm::vec3 bboxMin{bbox.first.x(), bbox.first.y(), bbox.first.z()};
    glm::vec3 bboxMax{bbox.second.x(), bbox.second.y(), bbox.second.z()};
    // register the grid
    polyscope::VolumeGrid* psGrid =
        polyscope::registerVolumeGrid("sample grid", {dimX, dimY, dimZ}, bboxMin, bboxMax);
    
    polyscope::VolumeGridNodeScalarQuantity* scalarQ = psGrid->addNodeScalarQuantity("node scalar1", g);


    scalarQ->setGridcubeVizEnabled(false);  // hide the default grid viz
    scalarQ->setIsosurfaceLevel(iso);        // set which isosurface we will visualize
    scalarQ->setIsosurfaceVizEnabled(true); // extracts the isosurface
    polyscope::show();

    // add a slice plane to cut through the grid while leaving the isosurface
    // untouched, as in the screenshot above
    scalarQ->setGridcubeVizEnabled(true);
    polyscope::SlicePlane* plane = polyscope::addSceneSlicePlane();
    scalarQ->setSlicePlanesAffectIsosurface(false);
    polyscope::show();

    // extract the isosurface as its own mesh structure
    scalarQ->registerIsosurfaceAsMesh("my isosurface mesh");


    delete cloudP;
    delete geoP;
}