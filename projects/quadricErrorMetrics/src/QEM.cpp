#include "QEM.h"
#include <queue>

 QEM::QEM(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo)
    : mesh(surfaceMesh), geometry(geo), errorMat(*surfaceMesh)
 {
     
     
 }

void QEM::update(size_t s) 
{
     
     size_t remain = mesh->nEdges()  - s;
     assert(remain >= 3); // must have a polygon

     
     using errPair = std::pair<double, Edge>;
     std::priority_queue<errPair, std::vector<errPair>, std::greater<errPair>> errors;

     EdgeData<bool> toReCalc(*mesh, false);

     // compute error matrix for every vertex
     for (Vertex v : mesh->vertices()) {
         errorMat[v] = ComputeErrorMat(v);
     }

     // select all valid pair (for simpilify, we just add all edges)
     // push all valid pair into the heap
     for (Edge e : mesh->edges()) {
         errors.push ({ ComputeLeastError(e).first, e});
     }

     // iteratively remove pair
     while (mesh->nEdges() > remain)
     {
         errPair err = errors.top();
         errors.pop();
         Edge e = err.second;
         if(e.isDead()) continue;
         if (toReCalc[e])
         {
             errors.push({ComputeLeastError(e).first, e});
             toReCalc[e] = false;
             continue;
         }

         // remove and updata position and error Matrix
         Vertex tail = e.firstVertex(), tip = e.secondVertex();
       
         Vector3 newPosition = ComputeLeastError(e).second; // use edge e before collapse

         Vertex RemainV = mesh->collapseEdgeTriangular(e);
         if (RemainV.getMesh()) {
             geometry->inputVertexPositions[RemainV] = newPosition;
             errorMat[RemainV] = ComputeErrorMat(RemainV);

             // recalulate error of vailid pairs of involved
             for (Edge e : RemainV.adjacentEdges()) {
                 toReCalc[e] = true;
             }
         }

     }

    
     
 }

 QEM::Mat QEM::ComputeErrorMat(Vertex v) 
 {
     Mat ret = Mat::Zero();
     Vector3 vposition = geometry->inputVertexPositions[v];
     for (Face f : v.adjacentFaces())
     {
         Vector3 fnormal = geometry->faceNormal(f);
         double d = -dot(vposition ,fnormal);
         Vec p{fnormal.x, fnormal.y, fnormal.z, d};
         Mat Kp = p * p.transpose();
         ret += Kp;
     }
     return ret;
 }

std::pair<double, Vector3> QEM::ComputeLeastError(Edge e) 
{
     Vector3 tip = geometry->inputVertexPositions[e.secondVertex()];
     Vector3 tail = geometry->inputVertexPositions[e.firstVertex()];
     Vector3 mid = (tip + tail)/2.;

     
     Vec midCoord{mid.x, mid.y, mid.z, 1};

     Mat quadricMat = errorMat[e.firstVertex()] + errorMat[e.secondVertex()];
     Mat auxQuad = quadricMat;
     quadricMat(3, 0) = 0.0f;
     quadricMat(3, 1) = 0.0f;
     quadricMat(3, 2) = 0.0f;
     quadricMat(3, 3) = 1.0f;
     Eigen::FullPivLU<Mat> lu(quadricMat);
     Vector3 ret;
     Vec b{0, 0, 0, 1.};
     double err = 0.;

     if (lu.isInvertible())
     {
         Vec slo = lu.solve(b);
         ret = Vector3{slo.x(), slo.y(), slo.z()};
         err = slo.transpose() * auxQuad * slo;
     } 
     else {
         ret = mid;
         err = midCoord.transpose() * auxQuad * midCoord;

     }

     return {err, ret};
 }

