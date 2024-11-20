#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <set>

struct Vertex {
    double x, y, z;
};

struct Face {
    int v1, v2, v3;
};

struct Edge {
    int v1, v2;

    bool operator<(const Edge& other) const {
        if (v1 == other.v1) return v2 < other.v2;
        return v1 < other.v1;
    }

    bool operator==(const Edge& other) const {
        return (v1 == other.v1 && v2 == other.v2) || (v1 == other.v2 && v2 == other.v1);
    }
};

// Хеш-функция для структуры Edge
struct EdgeHash {
    std::size_t operator()(const Edge& edge) const {
        return std::hash<int>()(edge.v1) ^ std::hash<int>()(edge.v2);
    }
};

void readObjFile(const std::string& filename, std::vector<Vertex>& vertices, std::vector<Face>& faces) {
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
                std::cout << ' ' << faces[i].v1 << ' ' << faces[i].v2 << ' ' << faces[i].v3 << std::endl;
            }
        }
    }

    file.close();
}


// std::vector<Edge> findEdges(const std::vector<Face>& faces) {
//     std::vector<Edge> edges;
//     std::unordered_map<Edge, int, EdgeHash> edgeCount;

//     for (const auto& face : faces) {
//         Edge e1 = {face.v1, face.v2};
//         Edge e2 = {face.v2, face.v3};
//         Edge e3 = {face.v3, face.v1};

//         edges.push_back(e1);
//         edges.push_back(e2);
//         edges.push_back(e3);

//         edgeCount[e1]++;
//         edgeCount[e2]++;
//         edgeCount[e3]++;
//     }

//     std::vector<Edge> sharedEdges;
//     for (const auto& edge : edges) {
//         if (edgeCount[edge] == 2) {
//             sharedEdges.push_back(edge);
//         }
//     }

//     return sharedEdges;
// }

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

void writeSharedEdgesToFile(const std::string& filename, const std::vector<Vertex>& vertices, const std::vector<Face>& faces, const std::unordered_map<Edge, std::vector<int>, EdgeHash>& edgeMap) {
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
                file << "Face " << faceIndex + 1 << ": " << face.v1 + 1 << " " << face.v2 + 1 << " " << face.v3 + 1 << "\n";
                file << "Vertices:\n";
                for (int vertexIndex : {face.v1, face.v2, face.v3}) {
                    const Vertex& vertex = vertices[vertexIndex];
                    file << std::setprecision(max_precision) << vertex.x << " " << vertex.y << " " << vertex.z << "\n";
                }
            }
            file << "\n";
        }
    }

    file.close();
}

int main() {
    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    readObjFile("angle_20x20.obj", vertices, faces);

    std::unordered_map<Edge, std::vector<int>, EdgeHash> edgeMap = findEdges(faces);

    std::cout << "Vertices:" << std::endl;
    int output_limit = 5;
    for (const auto& vertex : vertices) {
        if (output_limit <= 0) {
            break;
        }
        output_limit-=1;
        std::cout << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
    }

    output_limit = 5;
    std::cout << "Faces:" << std::endl;
    for (const auto& face : faces) {
        if (output_limit <= 0) {
            break;
        }
        output_limit-=1;
        std::cout << face.v1 << " " << face.v2 << " " << face.v3 << std::endl;
    }

    // std::cout << "Shared Edges:" << std::endl;
    // for (const auto& edge : sharedEdges) {
    //     std::cout << "Edge: " << edge.v1<< " " << edge.v2<< std::endl;
    // }

    // Нужен список смежных ребер для каждого треугольника
    output_limit = 5;
    std::cout << "Edges that belong to exactly two faces and their corresponding faces:" << std::endl;
    for (const auto& pair : edgeMap) {
        if (output_limit <= 0) {
            break;
        }
        output_limit-=1;
        const Edge& edge = pair.first;
        const std::vector<int>& faceIndices = pair.second;
        if (faceIndices.size() == 2) {
            std::cout << "Edge: " << edge.v1 + 1 << " " << edge.v2 + 1 << " belongs to faces: ";
            for (int faceIndex : faceIndices) {
                const Face& face = faces[faceIndex];
                std::cout << "Face " << faceIndex + 1 << ": " << face.v1 + 1 << " " << face.v2 + 1 << " " << face.v3 + 1 << "; ";
            }
            std::cout << std::endl;
        }
    }

    writeSharedEdgesToFile("shared_edges.txt", vertices, faces, edgeMap);

    return 0;
}