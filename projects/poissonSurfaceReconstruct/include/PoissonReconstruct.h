#include "geometrycentral/pointcloud/point_cloud.h"
#include "geometrycentral/pointcloud/point_cloud_io.h"
#include "geometrycentral/pointcloud/point_position_geometry.h"

using namespace geometrycentral;
using namespace geometrycentral::pointcloud;

class PoissonRec
{
public:
    // poisson reconstruction
    //
    // Inputs:
    //  P n by 3 points 
    //  N n by 3 normals of points
    //  dims 1 by 3 the dim of three direction
    //  bbox the bounding box of the grid
    // Outputs:
    //  x n by 3 the position of all grid positions
    //  g the scaler correspond to every grid nodes
  static void poisson_rec(const Eigen::MatrixXd& P, const Eigen::MatrixXd& N, Eigen::Vector3i& dims,
                          std::pair<Eigen::Vector3d, Eigen::Vector3d>& bbox,double& iso, Eigen::MatrixXd& x, Eigen::VectorXd& g);
    // Construct a gradient matrix for a finite-difference grid
    //
    // Inputs:
    //   nx  number of grid steps along the x-direction
    //   ny  number of grid steps along the y-direction
    //   nz  number of grid steps along the z-direction
    //   h  grid step size
    // Outputs:
    //   G  (nx-1)*ny*nz+ nx*(ny-1)*nz+ nx*ny*(nz-1) by nx*ny*nz sparse gradient
    //     matrix: G = [Dx;Dy;Dz]
    //
    static void fd_grad(const int nx, const int ny, const int nz, const double h, Eigen::SparseMatrix<double,Eigen::RowMajor>& G);
    // Construct a matrix of trilinear interpolation weights for a
    // finite-difference grid at a given set of points
    //
    // Inputs:
    //   nx  number of grid steps along the x-direction
    //   ny  number of grid steps along the y-direction
    //   nz  number of grid steps along the z-direction
    //   h  grid step size
    //   corner  list of bottom-left-front corner position of grid
    //   P  n by 3 list of query point locations
    // Outputs:
    //   W  n by (nx*ny*nz) sparse weights matrix
    static void fd_interpolate(const int nx, const int ny, const int nz, const double h, const Eigen::RowVector3d& corner,
                        const Eigen::MatrixXd& P, Eigen::SparseMatrix<double>& W);
    // Construct a partial derivative matrix for a finite-difference grid in a
    // given direction. Derivative are computed using first-order differences onto
    // a staggered grid
    //
    // Inputs:
    //   nx  number of grid steps along the x-direction
    //   ny  number of grid steps along the y-direction
    //   nz  number of grid steps along the z-direction
    //   h  grid step size
    //   dir  index indicating direction: 0-->x, 1-->y, 2-->z
    // Outputs:
    //   D  m by nx*ny*nz sparse partial derivative matrix, where:
    //     m = (nx-1)*ny*nz  if dir = 0
    //     m = nx*(ny-1)*nz  if dir = 1
    //     m = nx*ny*(nz-1)  otherwise (if dir = 2)
    //
    static void fd_partial_derivative(const int nx, const int ny, const int nz, const double h, const int dir,
                               Eigen::SparseMatrix<double>& D);
};