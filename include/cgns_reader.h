#pragma once

#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <array>

namespace cgns2openfoam {

struct Point3D {
    double x, y, z;
    
    bool operator==(const Point3D& other) const {
        const double tol = 1e-10;
        return std::abs(x - other.x) < tol && 
               std::abs(y - other.y) < tol && 
               std::abs(z - other.z) < tol;
    }
};

struct BoundaryPatch {
    std::string name;
    std::string type;  // wall, inlet, outlet, symmetry, cyclicAMI, etc.
    int startFace;
    int nFaces;
    
    // For cyclic/periodic patches
    std::string neighbourPatch;
    std::array<double, 3> rotationAxis = {0, 0, 0};
    std::array<double, 3> rotationCenter = {0, 0, 0};
    double rotationAngle = 0;  // in radians
    std::array<double, 3> separationVector = {0, 0, 0};
};

struct Zone {
    std::string name;
    std::vector<int> indices;
};

struct Face {
    std::vector<int> nodes;
    int owner = -1;
    int neighbour = -1;
    
    std::vector<int> sortedNodes() const {
        std::vector<int> sorted = nodes;
        std::sort(sorted.begin(), sorted.end());
        return sorted;
    }
};

// Interface connection between zones
struct ZoneInterface {
    std::string name;
    int zoneIndex1;
    int zoneIndex2;
    std::string zoneName1;
    std::string zoneName2;
    std::vector<int> faceIndices1;  // element indices in zone1
    std::vector<int> faceIndices2;  // element indices in zone2
    
    // Periodic/rotation info
    bool isPeriodic = false;
    std::array<double, 3> rotationAxis = {0, 0, 0};
    std::array<double, 3> rotationCenter = {0, 0, 0};
    double rotationAngle = 0;  // degrees
    std::array<double, 3> translation = {0, 0, 0};
};

struct MeshData {
    std::vector<Point3D> points;
    
    // All faces (internal first, then boundary)
    std::vector<Face> faces;
    int nInternalFaces = 0;
    
    // Boundary info
    std::vector<BoundaryPatch> boundaries;
    
    // Zones
    std::vector<Zone> cellZones;
    std::vector<Zone> faceZones;
    
    // Interfaces (for AMI/cyclic)
    std::vector<ZoneInterface> interfaces;
    
    int nCells = 0;
};

// Per-zone data during reading
struct ZoneData {
    std::string name;
    int cgnsIndex;
    
    std::vector<Point3D> points;
    std::vector<std::vector<int>> cellNodes;
    std::vector<int> cellTypes;
    
    // Face elements for BC
    std::map<int, std::vector<int>> elemIdxToFaceNodes;
    
    // BC info: elemIdx -> (bcName, bcType)
    std::map<int, std::pair<std::string, std::string>> bcInfo;
    
    // Offsets when merged
    int pointOffset = 0;
    int cellOffset = 0;
};

class CgnsReader {
public:
    bool open(const std::string& filepath);
    void close();

    bool readBase(int baseIndex, MeshData& mesh);
    bool isOpen() const { return fileId_ > 0; }
    int nBases() const { return nBases_; }

private:
    bool readZoneData(int baseIndex, int zoneIndex, ZoneData& zoneData);
    bool readZoneConnectivity(int baseIndex, int nZones, 
                              std::vector<ZoneData>& zones,
                              std::vector<ZoneInterface>& interfaces);
    void buildMultiZoneMesh(MeshData& mesh, 
                            std::vector<ZoneData>& zones,
                            const std::vector<ZoneInterface>& interfaces);
    
    int fileId_ = 0;
    int nBases_ = 0;
};

// Helper to get faces from cell type
std::vector<std::vector<int>> getCellFaces(int cellType, const std::vector<int>& nodes);

}  // namespace cgns2openfoam
