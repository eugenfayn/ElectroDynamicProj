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

class Triangle {
private:
    int a;  // first vertex of the shared edge
    int b;  // second vertex of the shared edge
    int c;  // vertex not on the shared edge
    std::array<std::array<double, 3>, 3> vertices;  // stores x,y,z coordinates for each vertex

public:
    Triangle();  // default constructor
    Triangle(int a_, int b_, int c_, 
            const double* v1, const double* v2, const double* v3);
    
    // Getters
    int getA() const;
    int getB() const;
    int getC() const;
    
    void getVertexCoordinates(int vertex, double& x, double& y, double& z) const;
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
