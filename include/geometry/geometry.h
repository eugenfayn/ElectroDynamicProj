#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <unordered_map>
#include <string>
#include <iomanip>
#include <iostream>
#include <limits>
#include <fstream>
#include <sstream>
#include <map>
#include <array>
#include <algorithm>
#include <complex>
#include <cmath>

class Vertex {
public:
    double x, y, z;
    int index;

    Vertex();
    Vertex(double x_, double y_, double z_, int index_);

    Vertex(const Vertex& other);
    
    bool operator==(const Vertex& other) const;
    bool operator<(const Vertex& other) const;
    Vertex operator+(const Vertex& other) const;
    Vertex operator-(const Vertex& other) const;
    Vertex operator/(double scalar) const;
    Vertex operator*(double scalar) const;
    Vertex operator*(const std::complex<double>& scalar) const;
    friend Vertex operator*(const std::complex<double>& scalar, const Vertex& vertex);
    Vertex& operator=(const Vertex& other);
    double scalar_product(const Vertex& other) const;
    double calc_S(const Vertex& b, const Vertex& c) const;
    double distance(const Vertex& other) const;
};

class Face {
public:
    int v1, v2, v3;

    Face();
    Face(int v1, int v2, int v3);
    bool operator==(const Face& other) const;
};

class Triangle {
private:
    Vertex a;  // first vertex of the shared edge
    Vertex b;  // second vertex of the shared edge
    Vertex c;  // vertex not on the shared edge
    int vertex_indices[3];

public:
    Triangle();  // default constructor
    Triangle(const Vertex& a_, const Vertex& b_, const Vertex& c_);
    
    // Getters
    const Vertex& getA() const { return a; }
    const Vertex& getB() const { return b; }
    const Vertex& getC() const { return c; }

    // Calculate area
    double calc_S() const;
    std::vector<int> getVertexIndices() const {
        // Return the vertex indices that make up this triangle
        return {this->a.index, this->b.index, this->c.index};
    }
};


// Function declaration for parsing shared edges
std::vector<Triangle> parseSharedEdges(const std::string& filename);

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

class MeshAnalyzer {
public:
    struct Edge {
        int v1, v2;
        
        Edge(int vertex1, int vertex2) {
            v1 = std::min(vertex1, vertex2);
            v2 = std::max(vertex1, vertex2);
        }
        
        bool operator<(const Edge& other) const {
            if (v1 != other.v1) return v1 < other.v1;
            return v2 < other.v2;
        }
    };

    static std::map<int, std::vector<int>> findNeighboringTriangles(const Triangle* triangles, const int N);
    static std::vector<Edge> getTriangleEdges(const Triangle& triangle);
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
