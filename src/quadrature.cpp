// quadrature.cpp
#include "quadrature/quadrature.h"
#include "math.h"

BarycentricCoord Quadrature::cartesianToBarycentric(
    const Vertex& point, 
    const Vertex& v1, 
    const Vertex& v2, 
    const Vertex& v3
) {
    double totalArea = v1.calc_S(v2, v3) * 2.0;
    
    BarycentricCoord result;
    result.xi1 = point.calc_S(v2, v3) * 2.0 / totalArea;
    result.xi2 = v1.calc_S(point, v3) * 2.0 / totalArea;
    result.xi3 = v1.calc_S(v2, point) * 2.0 / totalArea;
    
    return result;
}

std::vector<QuadraturePoint> Quadrature::getGauss3Points() {
    std::vector<QuadraturePoint> points(3);
    
    const double a = 1.0/6.0;
    const double b = 2.0/3.0;
    const double w = 1.0/3.0;
    
    points[0] = {b, a, a, w};
    points[1] = {a, b, a, w};
    points[2] = {a, a, b, w};
    
    return points;
}

std::vector<QuadraturePoint> Quadrature::getGauss4Points() {
    std::vector<QuadraturePoint> points(4);
    
    const double a = 0.585410196624969;
    const double b = 0.138196601125011;
    const double w = 0.25;
    
    points[0] = {a, b, b, w};
    points[1] = {b, a, b, w};
    points[2] = {b, b, a, w};
    points[3] = {1.0/3.0, 1.0/3.0, 1.0/3.0, w};
    
    return points;
}

Vertex Quadrature::barycentricToCartesian(
    const QuadraturePoint& bary,
    const Vertex& v1,
    const Vertex& v2,
    const Vertex& v3
) {
    return v1 * bary.xi1 + v2 * bary.xi2 + v3 * bary.xi3;
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
        // Use the function instead of inline calculation
        Vertex cartesianPoint = barycentricToCartesian(point, v1, v2, v3);
        result += func(cartesianPoint) * point.weight;
    }
    
    result *= area;
    return result;
}