#include "geometry/geometry.h"
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
        
        // Process edges
        auto edgeMap = findEdges(faces);
        
        // Print information
        printGeometryInfo(vertices, faces, edgeMap);
        
        // Write results
        writeSharedEdgesToFile(outputFile.string(), vertices, faces, edgeMap);
        
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
