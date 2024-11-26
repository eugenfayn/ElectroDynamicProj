#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <unordered_map>
#include <string>

class Vertex {
public:
    double x, y, z;

    Vertex();
    Vertex(double x, double y, double z);
    
    bool operator==(const Vertex& other) const;
    bool operator<(const Vertex& other) const;
    Vertex operator+(const Vertex& other) const;
    Vertex operator-(const Vertex& other) const;
    Vertex operator/(double scalar) const;
    Vertex operator*(double scalar) const;
    double scalar_product(const Vertex& other);
    double calc_S(const Vertex& b, const Vertex& c) const;
};

class Face {
public:
    int v1, v2, v3;

    Face();
    Face(int v1, int v2, int v3);
    bool operator==(const Face& other) const;
};

class Edge {
public:
    int v1, v2;

    Edge();
    Edge(int v1, int v2);
    bool operator<(const Edge& other) const;
    bool operator==(const Edge& other) const;
};

struct EdgeHash {
    std::size_t operator()(const Edge& edge) const;
};

// Geometry processing functions
void readGeometry(const std::string& filename, 
                 std::vector<Vertex>& vertices, 
                 std::vector<Face>& faces);

void printGeometryInfo(const std::vector<Vertex>& vertices, 
                      const std::vector<Face>& faces,
                      const std::unordered_map<Edge, std::vector<int>, EdgeHash>& edgeMap);

std::unordered_map<Edge, std::vector<int>, EdgeHash> findEdges(const std::vector<Face>& faces);

void writeSharedEdgesToFile(const std::string& filename, 
                          const std::vector<Vertex>& vertices, 
                          const std::vector<Face>& faces,
                          const std::unordered_map<Edge, std::vector<int>, EdgeHash>& edgeMap);

#endif
