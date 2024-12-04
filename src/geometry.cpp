#include "geometry/geometry.h"
#include "quadrature/quadrature.h"


// Vertex implementations
Vertex::Vertex() : x(0.0), y(0.0), z(0.0), index(-1) {}

Vertex::Vertex(double x_, double y_, double z_, int index_) 
    : x(x_), y(y_), z(z_), index(index_) {}

// Keep the existing operator implementations but update them to handle the index
bool Vertex::operator==(const Vertex& other) const {
    return x == other.x && y == other.y && z == other.z && index == other.index;
}

bool Vertex::operator<(const Vertex& other) const {
    if (x == other.x) {
        if (y == other.y) {
            if (z == other.z) {
                return index < other.index;
            }
            return z < other.z;
        }
        return y < other.y;
    }
    return x < other.x;
}

Vertex Vertex::operator+(const Vertex& other) const {
    return Vertex(x + other.x, y + other.y, z + other.z, index);
}

Vertex Vertex::operator-(const Vertex& other) const {
    return Vertex(x - other.x, y - other.y, z - other.z, index);
}

Vertex Vertex::operator/(double scalar) const {
    return Vertex(x/scalar, y/scalar, z/scalar, index);
}

Vertex Vertex::operator*(double scalar) const {
    return Vertex(x * scalar, y * scalar, z * scalar, index);
}

double Vertex::scalar_product(const Vertex& other) const {
    return x * other.x + y * other.y + z * other.z;
}

double Vertex::calc_S(const Vertex& b, const Vertex& c) const {
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

Triangle::Triangle() : a(), b(), c() {}

Triangle::Triangle(const Vertex& a_, const Vertex& b_, const Vertex& c_) 
    : a(a_), b(b_), c(c_) {}

std::vector<Triangle> parseSharedEdges(const std::string& filename) {
    std::vector<Triangle> triangles;
    std::ifstream file(filename);
    std::string line;
    
    while (std::getline(file, line)) {
        if (line.find("Edge:") != std::string::npos) {
            // Parse shared edge vertices indices
            int sharedV1Idx, sharedV2Idx;
            sscanf(line.c_str(), "Edge: %d %d", &sharedV1Idx, &sharedV2Idx);
            
            // Skip "Belongs to faces:" line
            std::getline(file, line);
            
            // Read first face
            std::getline(file, line);
            int faceNum1, f1v1, f1v2, f1v3;
            sscanf(line.c_str(), "Face %d: %d %d %d", &faceNum1, &f1v1, &f1v2, &f1v3);
            
            // Skip "Vertices:" line
            std::getline(file, line);
            
            // Read coordinates and create Vertex objects for first face
            std::map<int, Vertex> vertexMap1;
            for(int i = 0; i < 3; i++) {
                std::getline(file, line);
                double x, y, z;
                sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
                if (i == 0) vertexMap1[f1v1] = Vertex(x, y, z, f1v1);
                if (i == 1) vertexMap1[f1v2] = Vertex(x, y, z, f1v2);
                if (i == 2) vertexMap1[f1v3] = Vertex(x, y, z, f1v3);
            }
            
            // Read second face
            std::getline(file, line);
            int faceNum2, f2v1, f2v2, f2v3;
            sscanf(line.c_str(), "Face %d: %d %d %d", &faceNum2, &f2v1, &f2v2, &f2v3);
            
            // Skip "Vertices:" line
            std::getline(file, line);
            
            // Read coordinates and create Vertex objects for second face
            std::map<int, Vertex> vertexMap2;
            for(int i = 0; i < 3; i++) {
                std::getline(file, line);
                double x, y, z;
                sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
                if (i == 0) vertexMap2[f2v1] = Vertex(x, y, z, f2v1);
                if (i == 1) vertexMap2[f2v2] = Vertex(x, y, z, f2v2);
                if (i == 2) vertexMap2[f2v3] = Vertex(x, y, z, f2v3);
            }
            
            // Find third vertices (not in shared edge)
            int thirdV1Idx;
            if (f1v1 != sharedV1Idx && f1v1 != sharedV2Idx) thirdV1Idx = f1v1;
            else if (f1v2 != sharedV1Idx && f1v2 != sharedV2Idx) thirdV1Idx = f1v2;
            else thirdV1Idx = f1v3;
            
            int thirdV2Idx;
            if (f2v1 != sharedV1Idx && f2v1 != sharedV2Idx) thirdV2Idx = f2v1;
            else if (f2v2 != sharedV1Idx && f2v2 != sharedV2Idx) thirdV2Idx = f2v2;
            else thirdV2Idx = f2v3;
            
            // Create triangles with Vertex objects
            Triangle t1(vertexMap1[sharedV1Idx], 
                      vertexMap1[sharedV2Idx], 
                      vertexMap1[thirdV1Idx]);
            triangles.push_back(t1);
            
            Triangle t2(vertexMap2[sharedV1Idx], 
                      vertexMap2[sharedV2Idx], 
                      vertexMap2[thirdV2Idx]);
            triangles.push_back(t2);
        }
    }
    return triangles;
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
                    if (vertexIndex != edge.v1 && vertexIndex != edge.v2) {
                        const Vertex& vertex = vertices[vertexIndex];
                        file << std::setprecision(max_precision) 
                            << vertex.x << " " << vertex.y << " " << vertex.z << " Index: " << vertexIndex + 1 << "\n";
                    }
                    else {
                        const Vertex& vertex = vertices[vertexIndex];
                        file << std::setprecision(max_precision) 
                            << vertex.x << " " << vertex.y << " " << vertex.z << " Index: " << vertexIndex + 1 << "\n";
                    }
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
