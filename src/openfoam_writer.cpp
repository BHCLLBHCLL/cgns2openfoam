#include "openfoam_writer.h"
#include <fstream>
#include <iostream>
#include <filesystem>
#include <iomanip>
#include <cmath>

namespace fs = std::filesystem;
namespace cgns2openfoam {

static void writeFoamHeader(std::ostream& os, const std::string& className, 
                            const std::string& object) {
    os << "/*--------------------------------*- C++ -*----------------------------------*\\\n"
       << "| =========                 |                                                 |\n"
       << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
       << "|  \\\\    /   O peration     | Website:  https://openfoam.org                  |\n"
       << "|   \\\\  /    A nd           | Version:  converted by cgns2openfoam            |\n"
       << "|    \\\\/     M anipulation  |                                                 |\n"
       << "\\*---------------------------------------------------------------------------*/\n"
       << "FoamFile\n"
       << "{\n"
       << "    version     2.0;\n"
       << "    format      ascii;\n"
       << "    class       " << className << ";\n"
       << "    object      " << object << ";\n"
       << "}\n"
       << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n";
}

bool OpenFoamWriter::write(const std::string& caseDir, const MeshData& mesh) {
    return writePolyMesh(caseDir, mesh);
}

bool OpenFoamWriter::writePolyMesh(const std::string& caseDir, const MeshData& mesh) {
    fs::path polyMeshDir = fs::path(caseDir) / "constant" / "polyMesh";
    fs::create_directories(polyMeshDir);

    size_t nTotalFaces = mesh.faces.size();
    size_t nInternalFaces = static_cast<size_t>(mesh.nInternalFaces);

    // points
    {
        std::ofstream f(polyMeshDir / "points");
        if (!f) {
            std::cerr << "Cannot create points file" << std::endl;
            return false;
        }
        writeFoamHeader(f, "vectorField", "points");
        f << mesh.points.size() << "\n(\n";
        f << std::fixed << std::setprecision(16);
        for (const auto& p : mesh.points) {
            f << "(" << p.x << " " << p.y << " " << p.z << ")\n";
        }
        f << ")\n";
        f << "\n// ************************************************************************* //\n";
    }

    // faces
    {
        std::ofstream f(polyMeshDir / "faces");
        if (!f) {
            std::cerr << "Cannot create faces file" << std::endl;
            return false;
        }
        writeFoamHeader(f, "faceList", "faces");
        f << nTotalFaces << "\n(\n";
        for (const auto& face : mesh.faces) {
            f << face.nodes.size() << "(";
            for (size_t i = 0; i < face.nodes.size(); ++i) {
                f << face.nodes[i];
                if (i + 1 < face.nodes.size()) f << " ";
            }
            f << ")\n";
        }
        f << ")\n";
        f << "\n// ************************************************************************* //\n";
    }

    // owner
    {
        std::ofstream f(polyMeshDir / "owner");
        if (!f) return false;
        writeFoamHeader(f, "labelList", "owner");
        
        f << "\n";
        f << nTotalFaces << "\n(\n";
        for (const auto& face : mesh.faces) {
            f << face.owner << "\n";
        }
        f << ")\n";
        f << "\n// ************************************************************************* //\n";
    }

    // neighbour
    {
        std::ofstream f(polyMeshDir / "neighbour");
        if (!f) return false;
        writeFoamHeader(f, "labelList", "neighbour");
        f << "\n";
        f << nInternalFaces << "\n(\n";
        for (size_t i = 0; i < nInternalFaces; ++i) {
            f << mesh.faces[i].neighbour << "\n";
        }
        f << ")\n";
        f << "\n// ************************************************************************* //\n";
    }

    // boundary
    {
        std::ofstream f(polyMeshDir / "boundary");
        if (!f) return false;
        writeFoamHeader(f, "polyBoundaryMesh", "boundary");
        
        if (mesh.boundaries.empty()) {
            f << "0\n()\n";
        } else {
            f << mesh.boundaries.size() << "\n(\n";
            for (const auto& patch : mesh.boundaries) {
                f << "    " << patch.name << "\n";
                f << "    {\n";
                f << "        type            " << patch.type << ";\n";
                f << "        nFaces          " << patch.nFaces << ";\n";
                f << "        startFace       " << patch.startFace << ";\n";
                
                // Handle cyclic/AMI patches
                if (patch.type == "cyclicAMI" || patch.type == "cyclic") {
                    if (!patch.neighbourPatch.empty()) {
                        f << "        neighbourPatch  " << patch.neighbourPatch << ";\n";
                    }
                    
                    // Check for rotation
                    double axisMag = std::sqrt(
                        patch.rotationAxis[0]*patch.rotationAxis[0] +
                        patch.rotationAxis[1]*patch.rotationAxis[1] +
                        patch.rotationAxis[2]*patch.rotationAxis[2]);
                    
                    if (axisMag > 1e-10 && std::abs(patch.rotationAngle) > 1e-10) {
                        f << "        rotationAxis    (" 
                          << patch.rotationAxis[0] << " " 
                          << patch.rotationAxis[1] << " " 
                          << patch.rotationAxis[2] << ");\n";
                        f << "        rotationCentre  (" 
                          << patch.rotationCenter[0] << " " 
                          << patch.rotationCenter[1] << " " 
                          << patch.rotationCenter[2] << ");\n";
                        f << "        rotationAngle   " << patch.rotationAngle << ";\n";
                    }
                    
                    // Check for translation
                    double sepMag = std::sqrt(
                        patch.separationVector[0]*patch.separationVector[0] +
                        patch.separationVector[1]*patch.separationVector[1] +
                        patch.separationVector[2]*patch.separationVector[2]);
                    
                    if (sepMag > 1e-10) {
                        f << "        separationVector (" 
                          << patch.separationVector[0] << " " 
                          << patch.separationVector[1] << " " 
                          << patch.separationVector[2] << ");\n";
                    }
                }
                
                f << "    }\n";
            }
            f << ")\n";
        }
        f << "\n// ************************************************************************* //\n";
    }

    // cellZones - required for dynamic mesh
    {
        std::ofstream f(polyMeshDir / "cellZones");
        if (!f) return false;
        writeFoamHeader(f, "regIOobject", "cellZones");
        
        if (mesh.cellZones.empty()) {
            f << "0\n()\n";
        } else {
            f << mesh.cellZones.size() << "\n(\n";
            for (const auto& zone : mesh.cellZones) {
                f << zone.name << "\n";
                f << "{\n";
                f << "    type cellZone;\n";
                f << "    cellLabels List<label> " << zone.indices.size() << "(";
                for (size_t i = 0; i < zone.indices.size(); ++i) {
                    if (i % 10 == 0) f << "\n        ";
                    f << zone.indices[i];
                    if (i + 1 < zone.indices.size()) f << " ";
                }
                f << "\n    );\n";
                f << "}\n\n";
            }
            f << ")\n";
        }
        
        f << "\n// ************************************************************************* //\n";
    }

    // faceZones - write empty for now
    {
        std::ofstream f(polyMeshDir / "faceZones");
        if (!f) return false;
        writeFoamHeader(f, "regIOobject", "faceZones");
        
        f << "0\n()\n";
        
        f << "\n// ************************************************************************* //\n";
    }

    std::cout << "\nWrote polyMesh to " << polyMeshDir << "\n";
    std::cout << "  Points: " << mesh.points.size() << "\n";
    std::cout << "  Cells: " << mesh.nCells << "\n";
    std::cout << "  Total faces: " << nTotalFaces << "\n";
    std::cout << "  Internal faces: " << nInternalFaces << "\n";
    std::cout << "  Boundaries: " << mesh.boundaries.size() << "\n";
    for (const auto& patch : mesh.boundaries) {
        std::cout << "    - " << patch.name << " (" << patch.type << "): " 
                  << patch.nFaces << " faces, startFace=" << patch.startFace << "\n";
    }
    std::cout << "  Cell zones: " << mesh.cellZones.size() << "\n";
    for (const auto& zone : mesh.cellZones) {
        std::cout << "    - " << zone.name << ": " << zone.indices.size() << " cells\n";
    }

    return true;
}

}  // namespace cgns2openfoam
