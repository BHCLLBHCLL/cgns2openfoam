#pragma once

#include "cgns_reader.h"
#include <string>

namespace cgns2openfoam {

class OpenFoamWriter {
public:
    bool write(const std::string& caseDir, const MeshData& mesh);
    bool writePolyMesh(const std::string& caseDir, const MeshData& mesh);
};

}  // namespace cgns2openfoam
