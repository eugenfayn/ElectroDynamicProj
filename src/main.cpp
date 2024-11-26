#include "geometry/geometry.h"
#include "quadrature/quadrature.h"
#include "config.h"
#include <iostream>
#include <string>
#include <filesystem>
#include <chrono>

namespace fs = std::filesystem;

fs::path getInputPath(const std::string& filename) {
    return fs::path(PROJECT_INPUT_DIR) / filename;
}

fs::path getOutputPath(const std::string& filename) {
    return fs::path(PROJECT_OUTPUT_DIR) / filename;
}

void printBuildInfo() {
    std::cout << "Optimization level: -" << OPTIMIZATION_LEVEL << std::endl;
}

int main() {
    printBuildInfo();

    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    
    // Set the seed for the random number generator
    srand(time(NULL));

    auto start = std::chrono::high_resolution_clock::now();
    try {
        // Create output directory if it doesn't exist
        fs::create_directories(PROJECT_OUTPUT_DIR);

        // Get file paths
        fs::path inputFile = getInputPath("angle_20x20.obj");
        fs::path outputFile = getOutputPath("shared_edges.txt");

        // Check if input file exists
        if (!fs::exists(inputFile)) {
            std::cerr << "Input file not found: " << inputFile << std::endl;
            return 1;
        }

        // Read geometry
        readGeometry(inputFile.string(), vertices, faces);

        std::cout << "Number of vertices: " << vertices.size() << std::endl;
        std::cout << "Number of faces: " << faces.size() << std::endl;
        
        // Process edges
        auto edgeMap = findEdges(faces);
        
        // Print information
        printGeometryInfo(vertices, faces, edgeMap);
        
        // Write results
        writeSharedEdgesToFile(outputFile.string(), vertices, faces, edgeMap);

        // Example: integrate x*y over each face using Gauss-3 quadrature
        std::vector<QuadraturePoint> quadPoints = Quadrature::getGauss4Points();

        outputFile = getOutputPath("mesh_data.txt");
        writeGeometryAndQuadratureToFile(outputFile.string(), vertices, faces, quadPoints);

        // Print quadrature points to verify
        std::cout << "\nQuadrature points:" << std::endl;
        for (const auto& p : quadPoints) {
            std::cout << "Point: (" << p.xi1 << ", " << p.xi2 << ", " << p.xi3 
                    << "), weight: " << p.weight << std::endl;
        }
        
        // Define the function to integrate (example: f(x,y) = x*y)
        auto testFunction = [](const Vertex& v) {
            if (v.x * v.y == 0) {
                std::cout << "Warning: x or y is zero" << std::endl;
                throw std::runtime_error("x or y is zero");
            }
            return v.x * v.y;
        };

        double totalIntegral = 0.0;
        for (size_t i = 0; i < faces.size(); ++i) {
            const Face& face = faces[i];
            double faceIntegral = Quadrature::integrateOverFace(
                vertices, 
                face, 
                quadPoints, 
                testFunction
            );
            
            // Detailed output for first few faces
            if (i < 3) {
                const Vertex& v1 = vertices[face.v1];
                const Vertex& v2 = vertices[face.v2];
                const Vertex& v3 = vertices[face.v3];
                
                std::cout << "\nFace " << i << " details:" << std::endl;
                std::cout << "Vertices: " << std::endl;
                std::cout << "v1: (" << v1.x << ", " << v1.y << ", " << v1.z << ")" << std::endl;
                std::cout << "v2: (" << v2.x << ", " << v2.y << ", " << v2.z << ")" << std::endl;
                std::cout << "v3: (" << v3.x << ", " << v3.y << ", " << v3.z << ")" << std::endl;
                std::cout << "Integral: " << faceIntegral << std::endl;
            }
            
            totalIntegral += faceIntegral;
        }
        
        std::cout << "Total integral over mesh: " << totalIntegral << std::endl;
        
        std::cout << "Processing completed successfully." << std::endl;
        std::cout << "Output written to: " << outputFile << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() / 1000.0 << " seconds" << std::endl;
    
    return 0;
}
