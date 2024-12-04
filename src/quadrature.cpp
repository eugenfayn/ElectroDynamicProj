// quadrature.cpp
#include "quadrature/quadrature.h"


// quadrature.cpp
std::vector<QuadraturePoint> Quadrature::getGauss3Points() {
    std::vector<QuadraturePoint> points(3);
    
    const double a = 1.0/6.0;
    const double b = 2.0/3.0;
    const double w = 1.0/3.0;  // 1/6 * 2
    
    points[0] = {a, a, w};     // (1/6, 1/6)
    points[1] = {b, a, w};     // (2/3, 1/6)
    points[2] = {a, b, w};     // (1/6, 2/3)
    
    return points;
}

std::vector<QuadraturePoint> Quadrature::getGauss4Points() {
    std::vector<QuadraturePoint> points(4);
    
    const double a = 1.0/3.0;
    const double b = 3.0/5.0;
    const double c = 1.0/5.0;
    const double w = 25.0/48.0;
    const double w2 = -9.0/16.0;
    
    points[0] = {a, a, w2};    // center point
    points[1] = {b, c, w};     // outer points
    points[2] = {c, b, w};
    points[3] = {c, c, w};
    
    return points;
}

Vertex Quadrature::referenceToCartesian(
    const QuadraturePoint& point,
    const Vertex& v1,
    const Vertex& v2,
    const Vertex& v3
) {
    // Transform from reference coordinates (ξ,η) to physical coordinates
    // x = x1 + (x2-x1)ξ + (x3-x1)η
    // y = y1 + (y2-y1)ξ + (y3-y1)η
    // z = z1 + (z2-z1)ξ + (z3-z1)η
    
    return v1 + (v2 - v1) * point.xi + (v3 - v1) * point.eta;
}

double Quadrature::integrateOverFace(
    const std::vector<Vertex>& vertices,
    const Face& face,
    const std::vector<QuadraturePoint>& quadPoints,
    std::function<double(const Vertex&)> func
) {
    const Vertex& v1 = vertices[face.v1];
    const Vertex& v2 = vertices[face.v2];
    const Vertex& v3 = vertices[face.v3];
    
    // Calculate area using cross product for 3D triangle
    Vertex edge1 = v2 - v1; 
    Vertex edge2 = v3 - v1;
    
    double nx = edge1.y * edge2.z - edge1.z * edge2.y;
    double ny = edge1.z * edge2.x - edge1.x * edge2.z;
    double nz = edge1.x * edge2.y - edge1.y * edge2.x;
    
    double area = 0.5 * sqrt(nx*nx + ny*ny + nz*nz);
    
    double result = 0.0;
    
    for (const auto& point : quadPoints) {
        Vertex cartesianPoint = referenceToCartesian(point, v1, v2, v3);
        result += func(cartesianPoint) * point.weight;
    }
    
    result *= area;
    return result;
}
