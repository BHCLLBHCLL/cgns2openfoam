#include "cgns_reader.h"
#include "openfoam_writer.h"
#include <iostream>
#include <cstdlib>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: cgns2openfoam <input.cgns> <output_case_dir>\n"
                  << "  Converts CGNS mesh to OpenFOAM polyMesh format.\n"
                  << "  Supports multi-zone meshes with periodic/rotating interfaces.\n";
        return EXIT_FAILURE;
    }

    const std::string inputFile = argv[1];
    const std::string outputDir = argv[2];

    std::cout << "========================================\n";
    std::cout << "CGNS to OpenFOAM Converter\n";
    std::cout << "========================================\n\n";
    std::cout << "Input:  " << inputFile << "\n";
    std::cout << "Output: " << outputDir << "\n\n";

    cgns2openfoam::CgnsReader reader;
    if (!reader.open(inputFile)) {
        return EXIT_FAILURE;
    }

    cgns2openfoam::MeshData mesh;
    if (!reader.readBase(1, mesh)) {
        std::cerr << "Failed to read mesh from base 1\n";
        reader.close();
        return EXIT_FAILURE;
    }
    reader.close();

    std::cout << "\n========================================\n";
    std::cout << "Mesh Summary\n";
    std::cout << "========================================\n";
    std::cout << "  Points:         " << mesh.points.size() << "\n";
    std::cout << "  Cells:          " << mesh.nCells << "\n";
    std::cout << "  Total faces:    " << mesh.faces.size() << "\n";
    std::cout << "  Internal faces: " << mesh.nInternalFaces << "\n";
    std::cout << "  Boundary faces: " << (mesh.faces.size() - mesh.nInternalFaces) << "\n";
    std::cout << "\n  Boundaries (" << mesh.boundaries.size() << "):\n";
    for (const auto& bc : mesh.boundaries) {
        std::cout << "    - " << bc.name << " (" << bc.type << "): " 
                  << bc.nFaces << " faces\n";
    }
    std::cout << "\n  Cell zones (" << mesh.cellZones.size() << "):\n";
    for (const auto& zone : mesh.cellZones) {
        std::cout << "    - " << zone.name << ": " << zone.indices.size() << " cells\n";
    }
    if (!mesh.interfaces.empty()) {
        std::cout << "\n  Interfaces (" << mesh.interfaces.size() << "):\n";
        for (const auto& iface : mesh.interfaces) {
            std::cout << "    - " << iface.name << ": " << iface.zoneName1 
                      << " <-> " << iface.zoneName2;
            if (iface.isPeriodic) {
                std::cout << " [periodic]";
            }
            std::cout << "\n";
        }
    }

    cgns2openfoam::OpenFoamWriter writer;
    if (!writer.write(outputDir, mesh)) {
        std::cerr << "Failed to write OpenFOAM case\n";
        return EXIT_FAILURE;
    }

    std::cout << "\n========================================\n";
    std::cout << "Conversion completed successfully!\n";
    std::cout << "========================================\n";
    return EXIT_SUCCESS;
}
