#pragma once
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

class QEM
{
public:
    using Vec = Eigen::Vector4d;
    using Mat = Eigen::Matrix4d;
    ManifoldSurfaceMesh* mesh;
    VertexPositionGeometry* geometry;

    VertexData<Mat> errorMat;


    
    

    QEM(){};
    QEM(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo);
    
    void update(size_t s);
  

private:
    
    Mat ComputeErrorMat(Vertex v);
    std::pair<double, Vector3> ComputeLeastError(Edge e);
  
  
};