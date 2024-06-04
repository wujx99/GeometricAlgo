#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_cloud_io.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/volume_grid.h"
#include "RBF.h"

using namespace geometrycentral;
using namespace geometrycentral::pointcloud;

static std::pair<glm::vec3, glm::vec3> computeMinAndMax(const Vector<Vector3>& data) {
    Vector3 retMin{0, 0, 0}, retMax{0, 0, 0};
    for (size_t i = 0; i < data.size(); ++i) {
        retMin.x = retMin.x <= data[i].x ? retMin.x : data[i].x;
        retMin.y = retMin.y <= data[i].y ? retMin.y : data[i].y;
        retMin.z = retMin.z <= data[i].z ? retMin.z : data[i].z;
        retMax.x = retMax.x >= data[i].x ? retMax.x : data[i].x;
        retMax.y = retMax.y >= data[i].y ? retMax.y : data[i].y;
        retMax.z = retMax.z >= data[i].z ? retMax.z : data[i].z;
    }
    Vector3 padding = norm(retMin - retMax) * Vector3::constant(0.1);
    retMin -= padding;
    retMax += padding;
    return {{retMin.x, retMin.y, retMin.z}, {retMax.x, retMax.y, retMax.z}};
}

int main() {
    std::unique_ptr<PointCloud> cloud;
    std::unique_ptr<PointPositionGeometry> geom;
    std::string path = "D:/dev/GeometricAlgo/input/PointClouds/kitten_04.obj";
    std::tie(cloud, geom) = readPointCloud(path);

    PointCloud* cloudP = cloud.release();
    PointPositionGeometry* geoP = geom.release();

    auto minAndMax = computeMinAndMax(geoP->positions.raw());

    polyscope::init();

    // Register a point cloud
    // `points` is a Nx3 array-like container of points
    polyscope::registerPointCloud("my points", geoP->positions);
    RBFrecons rbf(cloudP, geoP);
    Vector3 p{0., 0., 1.};
    double ret = rbf.slove(p);

    uint32_t dimX = 100;
    uint32_t dimY = 100;
    uint32_t dimZ = 100;
    // register the grid
    polyscope::VolumeGrid* psGrid =
        polyscope::registerVolumeGrid("sample grid", {dimX, dimY, dimZ}, minAndMax.first, minAndMax.second);
    auto funcCall = [&](glm::vec3 position) {
        Vector3 p{position.x, position.y, position.z};
        return rbf.slove(p);
    };
    polyscope::VolumeGridNodeScalarQuantity* scalarQ = psGrid->addNodeScalarQuantityFromCallable("node scalar1", funcCall);
    
    
    scalarQ->setGridcubeVizEnabled(false);  // hide the default grid viz
    scalarQ->setIsosurfaceLevel(0.);       // set which isosurface we will visualize
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
