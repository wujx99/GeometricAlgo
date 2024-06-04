#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_cloud_io.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::pointcloud;

class RBFrecons
{
  public:
    using PointValPair = std::pair<Vector3,double>; 
    RBFrecons(PointCloud* c, PointPositionGeometry* g);
    double slove(Vector3 point);

    ~RBFrecons();
  private:
    Vector<PointValPair> allPointValPairs;
    

    Vector<PointValPair> coeff; // the coeff of all points
    double constCoeff{0.};
    Vector3 polyCoeff{0., 0., 0.};
    //Eigen::Vector4d polyCoeff; // the coeff of polynomial
    
    void addOffSurfacePoints(PointCloud* c, PointPositionGeometry* g, double dis = 1., size_t ratio = 2);
    void computeCoeff();
};