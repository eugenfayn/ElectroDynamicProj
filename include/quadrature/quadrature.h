// quadrature.h
#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <vector>
#include <functional>
#include "geometry/geometry.h"

struct BarycentricCoord {
    double xi1, xi2, xi3;
};

struct QuadraturePoint {
    double xi1, xi2, xi3;
    double weight;
};

class Quadrature {
public:
    static BarycentricCoord cartesianToBarycentric(
        const Vertex& point, 
        const Vertex& v1, 
        const Vertex& v2, 
        const Vertex& v3
    );

    static std::vector<QuadraturePoint> getGauss3Points();
    static std::vector<QuadraturePoint> getGauss4Points();

    static Vertex barycentricToCartesian(
        const QuadraturePoint& bary,
        const Vertex& v1,
        const Vertex& v2,
        const Vertex& v3
    );

    static double integrateOverFace(
        const std::vector<Vertex>& vertices,
        const Face& face,
        const std::vector<QuadraturePoint>& quadPoints,
        std::function<double(const Vertex&)> func
    );
};

void writeGeometryAndQuadratureToFile(const std::string& filename,
    const std::vector<Vertex>& vertices,
    const std::vector<Face>& faces,
    const std::vector<QuadraturePoint>& quadPoints);

#endif // QUADRATURE_H
