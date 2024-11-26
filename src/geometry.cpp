#include "geometry/geometry.h"
#include "quadrature/quadrature.h"
#include <iomanip>
#include <iostream>
#include <limits>
#include <fstream>
#include <sstream>

// Vertex implementations
Vertex::Vertex() : x(0.0), y(0.0), z(0.0) {}
Vertex::Vertex(double x, double y, double z) : x(x), y(y), z(z) {}

bool Vertex::operator==(const Vertex& other) const {
    return x == other.x && y == other.y && z == other.z;
}

bool Vertex::operator<(const Vertex& other) const {
    if (x == other.x) {
        if (y == other.y) {
            return z < other.z;
        }
        return y < other.y;
    }
    return x < other.x;
}

Vertex Vertex::operator+(const Vertex& other) const {
    return Vertex(x + other.x, y + other.y, z + other.z);
}

Vertex Vertex::operator-(const Vertex& other) const {
    return Vertex(x - other.x, y - other.y, z - other.z);
}

Vertex Vertex::operator/(double scalar) const {
    return Vertex(x/scalar, y/scalar, z/scalar);
}

Vertex Vertex::operator*(double scalar) const {
    return Vertex(x * scalar, y * scalar, z * scalar);
}

double Vertex::scalar_product(const Vertex& other) {
    return x * other.x + y * other.y + z * other.z;
}

double Vertex::calc_S(const Vertex& b, const Vertex& c) const{
    double S = 0.0;
    S = S + (b.x - x) * (c.y - y);
    S = S - (c.x - x) * (b.y - y);
    S = 0.5 * S;
    return S;
}

// Face implementations
Face::Face() : v1(0), v2(0), v3(0) {}
Face::Face(int v1, int v2, int v3) : v1(v1), v2(v2), v3(v3) {}

bool Face::operator==(const Face& other) const {
    return v1 == other.v1 && v2 == other.v2 && v3 == other.v3;
}

// Edge implementations
Edge::Edge() : v1(0), v2(0) {}
Edge::Edge(int v1, int v2) : v1(v1), v2(v2) {}

bool Edge::operator<(const Edge& other) const {
    if (v1 == other.v1) return v2 < other.v2;
    return v1 < other.v1;
}

bool Edge::operator==(const Edge& other) const {
    return (v1 == other.v1 && v2 == other.v2) || (v1 == other.v2 && v2 == other.v1);
}

std::size_t EdgeHash::operator()(const Edge& edge) const {
    return std::hash<int>()(edge.v1) ^ std::hash<int>()(edge.v2);
}

// Geometry processing functions implementations
void readGeometry(const std::string& filename, 
                 std::vector<Vertex>& vertices, 
                 std::vector<Face>& faces) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return;
    }

    std::string line;
    int numVertices = 0;
    int numFaces = 0;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "Points") {
            iss >> numVertices;
            vertices.resize(numVertices);
            for (int i = 0; i < numVertices; ++i) {
                std::getline(file, line);
                std::istringstream viss(line);
                viss >> vertices[i].x >> vertices[i].y >> vertices[i].z;
            }
        } else if (token == "Frames") {
            iss >> numFaces;
            faces.resize(numFaces);
            for (int i = 0; i < numFaces; ++i) {
                std::getline(file, line);
                std::istringstream fiss(line);
                fiss >> faces[i].v1 >> faces[i].v2 >> faces[i].v3;
                faces[i].v1 -= 1;
                faces[i].v2 -= 1;
                faces[i].v3 -= 1;
            }
        }
    }
    file.close();
}

void printGeometryInfo(const std::vector<Vertex>& vertices, 
                      const std::vector<Face>& faces,
                      const std::unordered_map<Edge, std::vector<int>, EdgeHash>& edgeMap) {
    // Print first 5 vertices
    std::cout << "Vertices:" << std::endl;
    for (size_t i = 0; i < std::min(vertices.size(), size_t(5)); ++i) {
        const auto& vertex = vertices[i];
        std::cout << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
    }

    // Print first 5 faces
    std::cout << "Faces:" << std::endl;
    for (size_t i = 0; i < std::min(faces.size(), size_t(5)); ++i) {
        const auto& face = faces[i];
        std::cout << face.v1 << " " << face.v2 << " " << face.v3 << std::endl;
    }

    // Print first 5 shared edges
    std::cout << "Edges that belong to exactly two faces and their corresponding faces:" << std::endl;
    size_t count = 0;
    for (const auto& [edge, faceIndices] : edgeMap) {
        if (count >= 5) break;
        if (faceIndices.size() == 2) {
            std::cout << "Edge: " << edge.v1 + 1 << " " << edge.v2 + 1 
                     << " belongs to faces: ";
            for (int faceIndex : faceIndices) {
                const Face& face = faces[faceIndex];
                std::cout << "Face " << faceIndex + 1 << ": " 
                         << face.v1 + 1 << " " << face.v2 + 1 << " " 
                         << face.v3 + 1 << "; ";
            }
            std::cout << std::endl;
            ++count;
        }
    }
}

std::unordered_map<Edge, std::vector<int>, EdgeHash> findEdges(const std::vector<Face>& faces) {
    std::unordered_map<Edge, std::vector<int>, EdgeHash> edgeMap;

    for (size_t i = 0; i < faces.size(); ++i) {
        const Face& face = faces[i];
        Edge e1 = {face.v1, face.v2};
        Edge e2 = {face.v2, face.v3};
        Edge e3 = {face.v3, face.v1};

        edgeMap[e1].push_back(i);
        edgeMap[e2].push_back(i);
        edgeMap[e3].push_back(i);
    }

    return edgeMap;
}

void writeSharedEdgesToFile(const std::string& filename, 
                          const std::vector<Vertex>& vertices, 
                          const std::vector<Face>& faces,
                          const std::unordered_map<Edge, std::vector<int>, EdgeHash>& edgeMap) {
    constexpr auto max_precision{std::numeric_limits<long double>::digits10 + 1};
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return;
    }

    for (const auto& pair : edgeMap) {
        const Edge& edge = pair.first;
        const std::vector<int>& faceIndices = pair.second;
        if (faceIndices.size() == 2) {
            file << "Edge: " << edge.v1 + 1 << " " << edge.v2 + 1 << "\n";
            file << "Belongs to faces:\n";
            for (int faceIndex : faceIndices) {
                const Face& face = faces[faceIndex];
                file << "Face " << faceIndex + 1 << ": " 
                     << face.v1 + 1 << " " << face.v2 + 1 << " " 
                     << face.v3 + 1 << "\n";
                file << "Vertices:\n";
                for (int vertexIndex : {face.v1, face.v2, face.v3}) {
                    const Vertex& vertex = vertices[vertexIndex];
                    file << std::setprecision(max_precision) 
                         << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
                }
            }
            file << "\n";
        }
    }

    file.close();
}

void writeGeometryAndQuadratureToFile(
    const std::string& filename,
    const std::vector<Vertex>& vertices,
    const std::vector<Face>& faces,
    const std::vector<QuadraturePoint>& quadPoints
) {
    std::ofstream file(filename);
    
    // Write vertices
    file << "Vertices\n";
    file << vertices.size() << "\n";
    for (const auto& v : vertices) {
        file << v.x << " " << v.y << " " << v.z << "\n";
    }
    
    // Write faces
    file << "Faces\n";
    file << faces.size() << "\n";
    for (const auto& f : faces) {
        file << f.v1 << " " << f.v2 << " " << f.v3 << "\n";
    }
    
    // Write quadrature points
    file << "QuadraturePoints\n";
    file << quadPoints.size() << "\n";
    for (const auto& p : quadPoints) {
        file << p.xi << " " << p.eta << " " << p.weight << "\n";
    }
    
    file.close();
}
