#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "QEM.h"

#include "args/args.hxx"

using namespace geometrycentral;
using namespace geometrycentral::surface;

int main(int argc, char** argv) {

    

    // Load a manifold surface mesh from file
    std::unique_ptr<ManifoldSurfaceMesh> mesh;
    std::unique_ptr<VertexPositionGeometry> geometry;
    std::string path = "D:/dev/GeometricAlgo/input/bunny.obj";

    std::tie(mesh, geometry) = readManifoldSurfaceMesh(path);

    polyscope::init();
    polyscope::registerSurfaceMesh("input mesh", geometry->inputVertexPositions, mesh->getFaceVertexList());

    size_t s = mesh->nEdges()-100;
    ManifoldSurfaceMesh* meshPtr = mesh.release();
    VertexPositionGeometry* geoPtr = geometry.release();
    QEM qemSimpilify(meshPtr, geoPtr);
    qemSimpilify.update(s);
    
    meshPtr->compress(); // be careful!!!
    
    auto* psMesh = polyscope::registerSurfaceMesh("subdiv mesh", geoPtr->inputVertexPositions, meshPtr->getFaceVertexList());


    polyscope::show();

    return EXIT_SUCCESS;
}