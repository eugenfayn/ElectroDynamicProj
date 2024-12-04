#include "geometry/geometry.h"
#include "quadrature/quadrature.h"
#include "integral.cpp"
#include "config.h"
#include <iostream>
#include <string>
#include <filesystem>
#include <chrono>
#include <complex>

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
    std::vector<Triangle> triangles;
    std::unordered_map<Edge, std::vector<int>, EdgeHash> edgeMap;
    
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
        edgeMap = findEdges(faces);
        
        // Print information
        printGeometryInfo(vertices, faces, edgeMap);
        
        // Write results
        writeSharedEdgesToFile(outputFile.string(), vertices, faces, edgeMap);

        std::vector<Triangle> triangles = parseSharedEdges(outputFile.string());
        std::cout << "Number of triangles: " << triangles.size() << std::endl;

        int ctr = 0;
        // To print triangle vertices:
        for(const Triangle& triangle : triangles) {
            std::cout << "Triangle vertices:" << std::endl;
            std::cout << "Vertex " << triangle.getA().index 
                    << " (shared edge): (" 
                    << triangle.getA().x << ", " 
                    << triangle.getA().y << ", " 
                    << triangle.getA().z << ")" << std::endl;
            std::cout << "Vertex " << triangle.getB().index 
                    << " (shared edge): (" 
                    << triangle.getB().x << ", " 
                    << triangle.getB().y << ", " 
                    << triangle.getB().z << ")" << std::endl;
            std::cout << "Vertex " << triangle.getC().index 
                    << " (not shared): (" 
                    << triangle.getC().x << ", " 
                    << triangle.getC().y << ", " 
                    << triangle.getC().z << ")" << std::endl;
            std::cout << "------------------------" << std::endl;
            ctr ++;
            if (ctr == 6) break;
        }
    
        // Example: integrate x*y over each face using Gauss-3 quadrature
        std::vector<QuadraturePoint> quadPoints = Quadrature::getGauss4Points();

        outputFile = getOutputPath("mesh_data.txt");
        writeGeometryAndQuadratureToFile(outputFile.string(), vertices, faces, quadPoints);

        // Print quadrature points to verify
        std::cout << "\nQuadrature points:" << std::endl;
        for (const auto& p : quadPoints) {
            std::cout << "Point: (" << p.xi << ", " << p.eta 
                    << "), weight: " << p.weight << std::endl;
        }
        
        // Example integration over a face
        auto testFunction = [](const Vertex& v) { return v.x * v.x + v.y * v.y; };

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

        std::cout << "Start BUILD MATRIX" << std::endl;
    
        // Calculate N (number of shared edges)
        int N = triangles.size() / 2;
        std::cout << "Matrix size: " << N << "x" << N << std::endl;
        std::cout << "Triangles size: " << triangles.size() << std::endl;
        
        // Allocate memory for the matrix
        std::complex<double>** A = new std::complex<double>*[N];
        for(int i = 0; i < N; i++) {
            A[i] = new std::complex<double>[N + 1];
        }

        // Build the matrix using triangles.data()
        const Triangle* trianglesPtr = triangles.data();
        BuildMatrix(A, N, trianglesPtr);

        // Optional: Print first few elements of matrix to verify
        std::cout << "First few elements of matrix:" << std::endl;
        for(int i = 0; i < std::min(3, N); i++) {
            for(int j = 0; j < std::min(3, N+1); j++) {
                std::cout << "A[" << i << "][" << j << "] = " 
                        << A[i][j].real() << " + " << A[i][j].imag() << "i\t";
            }
            std::cout << std::endl;
        }

        std::complex<double>* b = new std::complex<double>[N];
        // Initialize all elements to 0+0i (optional)
        for(int i = 0; i < N; i++) {
            b[i] = std::complex<double>(0.0, 0.0);
        }

        const Vertex polarization(0.0, 1.0, 0.0, -1);
        const Vertex tension(-1.0, 0.0, 0.0, -2);

        BuildRightPart(b, N, trianglesPtr, polarization, tension);

        SolveSLE(A, N, b);

        // Clean up
        for(int i = 0; i < N; i++) {
            delete[] A[i];
        }
        delete[] A;

        delete[] b;
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
