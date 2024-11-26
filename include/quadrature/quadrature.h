// quadrature.h
#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <vector>
#include <functional>
#include "geometry/geometry.h"

struct QuadraturePoint {
    double xi;     // ξ coordinate
    double eta;    // η coordinate
    double weight;
};

class Quadrature {
public:
    static std::vector<QuadraturePoint> getGauss3Points();
    static std::vector<QuadraturePoint> getGauss4Points();
    
    static Vertex referenceToCartesian(
        const QuadraturePoint& point,
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
