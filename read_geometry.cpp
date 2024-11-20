#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

struct Vertex {
    double x, y, z;
};

struct Face {
    int v1, v2, v3;
};

struct FrameGroup {
    int startFrame;
    int numFrames;
    int materialIndex;
};

void readObjFile(const std::string& filename, std::vector<Vertex>& vertices, std::vector<Face>& faces, std::vector<FrameGroup>& frameGroups) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return;
    }

    std::string line;
    int numVertices = 0;
    int numFaces = 0;
    int numFrameGroups = 0;

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
            }
        } else if (token == "FrameGroups") {
            iss >> numFrameGroups;
            frameGroups.resize(numFrameGroups);
            for (int i = 0; i < numFrameGroups; ++i) {
                std::getline(file, line);
                std::istringstream giss(line);
                giss >> frameGroups[i].startFrame >> frameGroups[i].numFrames >> frameGroups[i].materialIndex;
            }
        }
    }

    file.close();
}

int main() {
    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::vector<FrameGroup> frameGroups;

    readObjFile("ElectroDynamicProj/angle_20x20.obj", vertices, faces, frameGroups);

    std::cout << "Vertices:" << std::endl;

    int output_limit = 5;
    for (const auto& vertex : vertices) {
        if (output_limit <= 0) {
            std::cout << "output limit 1 " << output_limit << std::endl;
            break;
        }
        output_limit-=1;
        std::cout << vertex.x << " " << vertex.y << " " << vertex.z << std::endl;
    }

    output_limit = 5;
    std::cout << "Faces:" << std::endl;
    for (const auto& face : faces) {
        if (output_limit <= 0) {
            std::cout << "output limit 2 " << output_limit << std::endl;
            break;
        }
        output_limit-=1;
        std::cout << face.v1 << " " << face.v2 << " " << face.v3 << std::endl;
    }

    std::cout << "FrameGroups:" << std::endl;
    for (const auto& frameGroup : frameGroups) {
        std::cout << frameGroup.startFrame << " " << frameGroup.numFrames << " " << frameGroup.materialIndex << std::endl;
    }

    return 0;
}