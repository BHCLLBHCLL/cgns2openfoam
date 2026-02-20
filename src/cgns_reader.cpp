#include "cgns_reader.h"
#include <cgnslib.h>
#include <iostream>
#include <cstring>
#include <cmath>
#include <regex>

namespace cgns2openfoam {

// Sanitize name for OpenFOAM compatibility (no spaces, special chars)
static std::string sanitizeName(const std::string& name) {
    std::string result = name;
    // Replace spaces and special characters with underscores
    for (char& c : result) {
        if (c == ' ' || c == '/' || c == '\\' || c == ':' || c == ';' || 
            c == '(' || c == ')' || c == '[' || c == ']' || c == '{' || c == '}') {
            c = '_';
        }
    }
    // Remove leading/trailing underscores
    while (!result.empty() && result.front() == '_') result.erase(0, 1);
    while (!result.empty() && result.back() == '_') result.pop_back();
    // Replace multiple underscores with single
    std::regex multi_underscore("_+");
    result = std::regex_replace(result, multi_underscore, "_");
    return result;
}

bool CgnsReader::open(const std::string& filepath) {
    if (cg_open(filepath.c_str(), CG_MODE_READ, &fileId_) != CG_OK) {
        std::cerr << "Error opening CGNS file: " << filepath << std::endl;
        std::cerr << "  " << cg_get_error() << std::endl;
        return false;
    }
    if (cg_nbases(fileId_, &nBases_) != CG_OK) {
        std::cerr << "Error reading bases" << std::endl;
        close();
        return false;
    }
    return true;
}

void CgnsReader::close() {
    if (fileId_ > 0) {
        cg_close(fileId_);
        fileId_ = 0;
        nBases_ = 0;
    }
}

static bool isFaceElement(ElementType_t elemType) {
    return elemType == TRI_3 || elemType == TRI_6 ||
           elemType == QUAD_4 || elemType == QUAD_8 || elemType == QUAD_9;
}

static bool isCellElement(ElementType_t elemType) {
    return elemType == TETRA_4 || elemType == TETRA_10 ||
           elemType == PYRA_5 || elemType == PYRA_14 ||
           elemType == PENTA_6 || elemType == PENTA_15 || elemType == PENTA_18 ||
           elemType == HEXA_8 || elemType == HEXA_20 || elemType == HEXA_27;
}

static std::string bcTypeToFoamType(BCType_t bcType) {
    switch (bcType) {
        case BCWall:
        case BCWallViscous:
        case BCWallViscousHeatFlux:
        case BCWallViscousIsothermal:
        case BCWallInviscid:
            return "wall";
        case BCInflow:
        case BCInflowSubsonic:
        case BCInflowSupersonic:
            return "inlet";
        case BCOutflow:
        case BCOutflowSubsonic:
        case BCOutflowSupersonic:
        case BCExtrapolate:
            return "outlet";
        case BCSymmetryPlane:
            return "symmetryPlane";
        case BCSymmetryPolar:
            return "symmetry";
        case BCFarfield:
            return "freestream";
        default:
            return "patch";
    }
}

std::vector<std::vector<int>> getCellFaces(int cellType, const std::vector<int>& nodes) {
    std::vector<std::vector<int>> faces;
    
    switch (cellType) {
        case TETRA_4:
        case TETRA_10:
            faces.push_back({nodes[0], nodes[2], nodes[1]});
            faces.push_back({nodes[0], nodes[1], nodes[3]});
            faces.push_back({nodes[1], nodes[2], nodes[3]});
            faces.push_back({nodes[0], nodes[3], nodes[2]});
            break;
            
        case PYRA_5:
        case PYRA_14:
            faces.push_back({nodes[0], nodes[3], nodes[2], nodes[1]});
            faces.push_back({nodes[0], nodes[1], nodes[4]});
            faces.push_back({nodes[1], nodes[2], nodes[4]});
            faces.push_back({nodes[2], nodes[3], nodes[4]});
            faces.push_back({nodes[3], nodes[0], nodes[4]});
            break;
            
        case PENTA_6:
        case PENTA_15:
        case PENTA_18:
            faces.push_back({nodes[0], nodes[2], nodes[1]});
            faces.push_back({nodes[3], nodes[4], nodes[5]});
            faces.push_back({nodes[0], nodes[1], nodes[4], nodes[3]});
            faces.push_back({nodes[1], nodes[2], nodes[5], nodes[4]});
            faces.push_back({nodes[2], nodes[0], nodes[3], nodes[5]});
            break;
            
        case HEXA_8:
        case HEXA_20:
        case HEXA_27:
            faces.push_back({nodes[0], nodes[3], nodes[2], nodes[1]});
            faces.push_back({nodes[4], nodes[5], nodes[6], nodes[7]});
            faces.push_back({nodes[0], nodes[1], nodes[5], nodes[4]});
            faces.push_back({nodes[2], nodes[3], nodes[7], nodes[6]});
            faces.push_back({nodes[1], nodes[2], nodes[6], nodes[5]});
            faces.push_back({nodes[0], nodes[4], nodes[7], nodes[3]});
            break;
            
        default:
            break;
    }
    
    return faces;
}

static std::string faceKey(std::vector<int> nodes) {
    std::sort(nodes.begin(), nodes.end());
    std::string key;
    for (int n : nodes) {
        key += std::to_string(n) + "_";
    }
    return key;
}

bool CgnsReader::readBase(int baseIndex, MeshData& mesh) {
    char basename[33];
    int cellDim, physDim;
    if (cg_base_read(fileId_, baseIndex, basename, &cellDim, &physDim) != CG_OK) {
        std::cerr << "Error reading base: " << cg_get_error() << std::endl;
        return false;
    }

    std::cout << "Reading base: " << basename << " (cellDim=" << cellDim << ", physDim=" << physDim << ")\n";

    int nZones;
    if (cg_nzones(fileId_, baseIndex, &nZones) != CG_OK) {
        std::cerr << "Error reading zones" << std::endl;
        return false;
    }
    
    std::cout << "  Found " << nZones << " zones\n";

    // Read all zones
    std::vector<ZoneData> zones(nZones);
    for (int iz = 1; iz <= nZones; ++iz) {
        if (!readZoneData(baseIndex, iz, zones[iz-1])) {
            return false;
        }
        zones[iz-1].cgnsIndex = iz;
    }

    // Read zone connectivity (interfaces)
    std::vector<ZoneInterface> interfaces;
    readZoneConnectivity(baseIndex, nZones, zones, interfaces);

    // Build merged mesh
    mesh.points.clear();
    mesh.faces.clear();
    mesh.boundaries.clear();
    mesh.cellZones.clear();
    mesh.faceZones.clear();
    mesh.interfaces.clear();
    mesh.nCells = 0;
    mesh.nInternalFaces = 0;

    buildMultiZoneMesh(mesh, zones, interfaces);
    
    mesh.interfaces = interfaces;

    return true;
}

bool CgnsReader::readZoneData(int baseIndex, int zoneIndex, ZoneData& zoneData) {
    char zoneName[33];
    cgsize_t size[9];
    ZoneType_t zoneType;
    
    if (cg_zone_type(fileId_, baseIndex, zoneIndex, &zoneType) != CG_OK) {
        std::cerr << "Error reading zone type" << std::endl;
        return false;
    }
    
    if (zoneType != Unstructured) {
        std::cerr << "Only unstructured zones are supported" << std::endl;
        return false;
    }
    
    if (cg_zone_read(fileId_, baseIndex, zoneIndex, zoneName, size) != CG_OK) {
        std::cerr << "Error reading zone " << zoneIndex << std::endl;
        return false;
    }

    zoneData.name = sanitizeName(zoneName);
    cgsize_t nVertex = size[0];
    cgsize_t nCell = size[1];
    
    std::cout << "  Zone " << zoneIndex << ": " << zoneName 
              << " (vertices=" << nVertex << ", cells=" << nCell << ")\n";

    // Read coordinates
    cgsize_t rmin = 1, rmax = nVertex;
    std::vector<double> x(nVertex), y(nVertex), z(nVertex);
    
    if (cg_coord_read(fileId_, baseIndex, zoneIndex, "CoordinateX", RealDouble, &rmin, &rmax, x.data()) != CG_OK) {
        std::cerr << "Error reading X coordinates" << std::endl;
        return false;
    }
    if (cg_coord_read(fileId_, baseIndex, zoneIndex, "CoordinateY", RealDouble, &rmin, &rmax, y.data()) != CG_OK) {
        std::cerr << "Error reading Y coordinates" << std::endl;
        return false;
    }
    if (cg_coord_read(fileId_, baseIndex, zoneIndex, "CoordinateZ", RealDouble, &rmin, &rmax, z.data()) != CG_OK) {
        std::cerr << "Error reading Z coordinates" << std::endl;
        return false;
    }

    zoneData.points.resize(nVertex);
    for (cgsize_t i = 0; i < nVertex; ++i) {
        zoneData.points[i] = {x[i], y[i], z[i]};
    }

    // Read element sections
    int nSections;
    if (cg_nsections(fileId_, baseIndex, zoneIndex, &nSections) != CG_OK) {
        std::cerr << "Error reading sections count" << std::endl;
        return false;
    }

    for (int is = 1; is <= nSections; ++is) {
        char sectName[33];
        ElementType_t elemType;
        cgsize_t start, end;
        int nBndry, parentFlag;
        
        if (cg_section_read(fileId_, baseIndex, zoneIndex, is,
                sectName, &elemType, &start, &end, &nBndry, &parentFlag) != CG_OK) {
            continue;
        }

        std::cout << "    Section " << is << ": " << sectName 
                  << " (type=" << elemType << ", range=" << start << "-" << end << ")\n";

        cgsize_t nElemInSect = end - start + 1;
        cgsize_t elemDataSize;
        if (cg_ElementDataSize(fileId_, baseIndex, zoneIndex, is, &elemDataSize) != CG_OK) {
            continue;
        }

        std::vector<cgsize_t> elemData(elemDataSize);
        if (cg_elements_read(fileId_, baseIndex, zoneIndex, is, elemData.data(), nullptr) != CG_OK) {
            continue;
        }

        int nNodesPerElem = 0;
        if (cg_npe(elemType, &nNodesPerElem) != CG_OK || nNodesPerElem <= 0) {
            continue;
        }

        if (isCellElement(elemType)) {
            size_t idx = 0;
            for (cgsize_t ie = 0; ie < nElemInSect && idx < elemData.size(); ++ie) {
                std::vector<int> nodes;
                for (int in = 0; in < nNodesPerElem && idx < elemData.size(); ++in) {
                    nodes.push_back(static_cast<int>(elemData[idx++] - 1));
                }
                zoneData.cellNodes.push_back(nodes);
                zoneData.cellTypes.push_back(elemType);
            }
        } else if (isFaceElement(elemType)) {
            size_t idx = 0;
            for (cgsize_t ie = 0; ie < nElemInSect && idx < elemData.size(); ++ie) {
                std::vector<int> nodes;
                for (int in = 0; in < nNodesPerElem && idx < elemData.size(); ++in) {
                    nodes.push_back(static_cast<int>(elemData[idx++] - 1));
                }
                int elemIdx = static_cast<int>(start + ie);
                zoneData.elemIdxToFaceNodes[elemIdx] = nodes;
                // Default: use section name as BC name
                zoneData.bcInfo[elemIdx] = {sanitizeName(sectName), "patch"};
            }
        }
    }

    // Read boundary conditions - these override section names
    int nBocos;
    if (cg_nbocos(fileId_, baseIndex, zoneIndex, &nBocos) != CG_OK) {
        nBocos = 0;
    }

    std::cout << "    Found " << nBocos << " boundary conditions\n";

    for (int ib = 1; ib <= nBocos; ++ib) {
        char bocoName[33];
        BCType_t bocoType;
        PointSetType_t ptsetType;
        cgsize_t npnts;
        int normalIndex[3];
        cgsize_t normalListSize;
        DataType_t normalDataType;
        int nDataSet;

        if (cg_boco_info(fileId_, baseIndex, zoneIndex, ib,
                bocoName, &bocoType, &ptsetType, &npnts,
                normalIndex, &normalListSize, &normalDataType, &nDataSet) != CG_OK) {
            continue;
        }

        std::cout << "      BC " << ib << ": " << bocoName << " (type=" << bocoType 
                  << ", ptset=" << ptsetType << ", npnts=" << npnts << ")\n";

        std::vector<cgsize_t> pnts(npnts);
        if (cg_boco_read(fileId_, baseIndex, zoneIndex, ib, pnts.data(), nullptr) != CG_OK) {
            continue;
        }

        std::string foamType = bcTypeToFoamType(bocoType);

        std::vector<int> bcElemIndices;
        if (ptsetType == PointRange) {
            for (cgsize_t ie = pnts[0]; ie <= pnts[1]; ++ie) {
                bcElemIndices.push_back(static_cast<int>(ie));
            }
        } else if (ptsetType == PointList || ptsetType == ElementList) {
            for (cgsize_t i = 0; i < npnts; ++i) {
                bcElemIndices.push_back(static_cast<int>(pnts[i]));
            }
        }
        
        for (int elemIdx : bcElemIndices) {
            zoneData.bcInfo[elemIdx] = {sanitizeName(bocoName), foamType};
        }
    }

    return true;
}

bool CgnsReader::readZoneConnectivity(int baseIndex, int nZones, 
                                      std::vector<ZoneData>& zones,
                                      std::vector<ZoneInterface>& interfaces) {
    for (int iz = 1; iz <= nZones; ++iz) {
        // Read 1-to-1 connectivity
        int n1to1;
        if (cg_n1to1(fileId_, baseIndex, iz, &n1to1) == CG_OK && n1to1 > 0) {
            std::cout << "    Zone " << iz << " has " << n1to1 << " 1-to-1 connections\n";
            
            for (int ic = 1; ic <= n1to1; ++ic) {
                char connName[33], donorName[33];
                cgsize_t range[6], donorRange[6];
                int transform[3];
                
                if (cg_1to1_read(fileId_, baseIndex, iz, ic, connName, donorName,
                        range, donorRange, transform) == CG_OK) {
                    std::cout << "      Conn " << ic << ": " << connName 
                              << " -> " << donorName << "\n";
                    
                    int donorZoneIdx = -1;
                    for (int dz = 0; dz < nZones; ++dz) {
                        if (zones[dz].name == donorName) {
                            donorZoneIdx = dz;
                            break;
                        }
                    }
                    
                    if (donorZoneIdx >= 0 && donorZoneIdx != iz - 1) {
                        ZoneInterface iface;
                        iface.name = sanitizeName(connName);
                        iface.zoneIndex1 = iz - 1;
                        iface.zoneIndex2 = donorZoneIdx;
                        iface.zoneName1 = zones[iz-1].name;
                        iface.zoneName2 = sanitizeName(donorName);
                        interfaces.push_back(iface);
                    }
                }
            }
        }
        
        // Read general connectivity
        int nConns;
        if (cg_nconns(fileId_, baseIndex, iz, &nConns) == CG_OK && nConns > 0) {
            std::cout << "    Zone " << iz << " has " << nConns << " general connections\n";
            
            for (int ic = 1; ic <= nConns; ++ic) {
                char connName[33], donorName[33];
                GridLocation_t location;
                GridConnectivityType_t connType;
                PointSetType_t ptsetType, donorPtsetType;
                cgsize_t npnts, donorDataSize;
                ZoneType_t donorZoneType;
                DataType_t donorDataType;
                
                if (cg_conn_info(fileId_, baseIndex, iz, ic, connName, &location,
                        &connType, &ptsetType, &npnts, donorName, &donorZoneType,
                        &donorPtsetType, &donorDataType, &donorDataSize) == CG_OK) {
                    
                    std::cout << "      Conn " << ic << ": " << connName 
                              << " -> " << donorName << " (type=" << connType << ")\n";
                    
                    int donorZoneIdx = -1;
                    for (int dz = 0; dz < nZones; ++dz) {
                        if (zones[dz].name == donorName) {
                            donorZoneIdx = dz;
                            break;
                        }
                    }
                    
                    if (donorZoneIdx < 0) continue;
                    
                    std::vector<cgsize_t> pnts(npnts);
                    std::vector<cgsize_t> donorPnts(donorDataSize);
                    
                    if (cg_conn_read(fileId_, baseIndex, iz, ic, pnts.data(),
                            donorDataType, donorPnts.data()) != CG_OK) {
                        continue;
                    }
                    
                    ZoneInterface iface;
                    iface.name = sanitizeName(connName);
                    iface.zoneIndex1 = iz - 1;
                    iface.zoneIndex2 = donorZoneIdx;
                    iface.zoneName1 = zones[iz-1].name;
                    iface.zoneName2 = sanitizeName(donorName);
                    
                    for (cgsize_t i = 0; i < npnts; ++i) {
                        iface.faceIndices1.push_back(static_cast<int>(pnts[i]));
                    }
                    for (cgsize_t i = 0; i < donorDataSize; ++i) {
                        iface.faceIndices2.push_back(static_cast<int>(donorPnts[i]));
                    }
                    
                    // Check for periodic
                    float rotCenter[3], rotAngle[3], translation[3];
                    if (cg_conn_periodic_read(fileId_, baseIndex, iz, ic, 
                            rotCenter, rotAngle, translation) == CG_OK) {
                        iface.isPeriodic = true;
                        iface.rotationCenter = {rotCenter[0], rotCenter[1], rotCenter[2]};
                        
                        double ax = rotAngle[0], ay = rotAngle[1], az = rotAngle[2];
                        double mag = std::sqrt(ax*ax + ay*ay + az*az);
                        if (mag > 1e-10) {
                            iface.rotationAxis = {ax/mag, ay/mag, az/mag};
                            iface.rotationAngle = mag * 180.0 / M_PI;
                        }
                        iface.translation = {translation[0], translation[1], translation[2]};
                        
                        std::cout << "        Periodic: rot=(" << rotAngle[0] << "," 
                                  << rotAngle[1] << "," << rotAngle[2] << ") trans=("
                                  << translation[0] << "," << translation[1] << "," 
                                  << translation[2] << ")\n";
                    }
                    
                    interfaces.push_back(iface);
                }
            }
        }
    }
    
    return true;
}

struct BCInfo {
    std::string name;
    std::string type;
};

void CgnsReader::buildMultiZoneMesh(MeshData& mesh, 
                                    std::vector<ZoneData>& zones,
                                    const std::vector<ZoneInterface>& interfaces) {
    // Calculate offsets and merge points
    int pointOffset = 0;
    int cellOffset = 0;
    
    for (auto& zone : zones) {
        zone.pointOffset = pointOffset;
        zone.cellOffset = cellOffset;
        
        for (const auto& p : zone.points) {
            mesh.points.push_back(p);
        }
        
        pointOffset += static_cast<int>(zone.points.size());
        cellOffset += static_cast<int>(zone.cellNodes.size());
    }
    
    mesh.nCells = cellOffset;
    
    // Build faceKey to BC info mapping (with global point indices)
    // Also track faceZone info
    std::unordered_map<std::string, BCInfo> faceKeyToBc;
    std::unordered_map<std::string, std::string> faceKeyToZoneName;  // for faceZones
    
    // Track interface faces with pairing info
    // key -> (patchName, neighbourPatchName, isPeriodic)
    struct InterfaceData {
        std::string patchName;
        std::string neighbourPatch;
        bool isPeriodic;
        std::array<double, 3> rotationAxis;
        std::array<double, 3> rotationCenter;
        double rotationAngle;
    };
    std::set<std::string> interfaceFaceKeys;
    std::map<std::string, InterfaceData> interfaceInfo;
    
    for (const auto& iface : interfaces) {
        auto& zone1 = zones[iface.zoneIndex1];
        auto& zone2 = zones[iface.zoneIndex2];
        
        std::string patchName1 = iface.name + "_" + zone1.name;
        std::string patchName2 = iface.name + "_" + zone2.name;
        
        for (int elemIdx : iface.faceIndices1) {
            auto it = zone1.elemIdxToFaceNodes.find(elemIdx);
            if (it != zone1.elemIdxToFaceNodes.end()) {
                std::vector<int> globalNodes;
                for (int n : it->second) {
                    globalNodes.push_back(n + zone1.pointOffset);
                }
                std::string key = faceKey(globalNodes);
                interfaceFaceKeys.insert(key);
                interfaceInfo[key] = {patchName1, patchName2, iface.isPeriodic,
                                      iface.rotationAxis, iface.rotationCenter, iface.rotationAngle};
            }
        }
        
        for (int elemIdx : iface.faceIndices2) {
            auto it = zone2.elemIdxToFaceNodes.find(elemIdx);
            if (it != zone2.elemIdxToFaceNodes.end()) {
                std::vector<int> globalNodes;
                for (int n : it->second) {
                    globalNodes.push_back(n + zone2.pointOffset);
                }
                std::string key = faceKey(globalNodes);
                interfaceFaceKeys.insert(key);
                // Reverse rotation for the other side
                interfaceInfo[key] = {patchName2, patchName1, iface.isPeriodic,
                                      iface.rotationAxis, iface.rotationCenter, -iface.rotationAngle};
            }
        }
    }
    
    // Build BC mapping for each zone (multi-zone: prepend zone name if needed)
    bool multiZone = zones.size() > 1;
    
    for (auto& zone : zones) {
        for (const auto& [elemIdx, nodes] : zone.elemIdxToFaceNodes) {
            std::vector<int> globalNodes;
            for (int n : nodes) {
                globalNodes.push_back(n + zone.pointOffset);
            }
            std::string key = faceKey(globalNodes);
            
            // Skip interface faces for BC assignment
            if (interfaceFaceKeys.count(key)) continue;
            
            auto bcIt = zone.bcInfo.find(elemIdx);
            if (bcIt != zone.bcInfo.end()) {
                std::string bcName = bcIt->second.first;
                // For multi-zone, prepend zone name to avoid name collisions
                if (multiZone) {
                    bcName = zone.name + "_" + bcName;
                }
                faceKeyToBc[key] = {bcName, bcIt->second.second};
                faceKeyToZoneName[key] = bcName;  // Track for faceZone
            }
        }
    }
    
    // Extract faces from all cells
    struct FaceInfo {
        std::vector<int> nodes;
        int owner = -1;
        int neighbour = -1;
        std::string bcName;
        std::string bcType;
        bool isInterface = false;
        std::string neighbourPatch;
        std::array<double, 3> rotationAxis = {0, 0, 0};
        std::array<double, 3> rotationCenter = {0, 0, 0};
        double rotationAngle = 0;
    };
    std::unordered_map<std::string, FaceInfo> faceMap;
    
    for (size_t zi = 0; zi < zones.size(); ++zi) {
        auto& zone = zones[zi];
        int nCells = static_cast<int>(zone.cellNodes.size());
        
        for (int ci = 0; ci < nCells; ++ci) {
            int globalCellIdx = ci + zone.cellOffset;
            
            std::vector<int> globalNodes;
            for (int n : zone.cellNodes[ci]) {
                globalNodes.push_back(n + zone.pointOffset);
            }
            
            auto faces = getCellFaces(zone.cellTypes[ci], globalNodes);
            
            for (auto& faceNodes : faces) {
                std::string key = faceKey(faceNodes);
                
                auto it = faceMap.find(key);
                if (it == faceMap.end()) {
                    FaceInfo info;
                    info.nodes = faceNodes;
                    info.owner = globalCellIdx;
                    faceMap[key] = info;
                } else {
                    if (globalCellIdx < it->second.owner) {
                        it->second.neighbour = it->second.owner;
                        it->second.owner = globalCellIdx;
                        it->second.nodes = faceNodes;
                    } else {
                        it->second.neighbour = globalCellIdx;
                    }
                }
            }
        }
    }
    
    // Apply BC and interface info
    for (auto& [key, info] : faceMap) {
        if (info.neighbour < 0) {
            auto bcIt = faceKeyToBc.find(key);
            if (bcIt != faceKeyToBc.end()) {
                info.bcName = bcIt->second.name;
                info.bcType = bcIt->second.type;
            }
            
            auto ifIt = interfaceInfo.find(key);
            if (ifIt != interfaceInfo.end()) {
                info.bcName = ifIt->second.patchName;
                info.bcType = "cyclicAMI";  // All interfaces use cyclicAMI
                info.isInterface = true;
                info.neighbourPatch = ifIt->second.neighbourPatch;
                info.rotationAxis = ifIt->second.rotationAxis;
                info.rotationCenter = ifIt->second.rotationCenter;
                info.rotationAngle = ifIt->second.rotationAngle;
            }
        }
    }
    
    // Separate internal and boundary faces
    std::vector<Face> internalFaces;
    
    // Store boundary info with cyclicAMI data
    struct BoundaryEntry {
        std::string bcType;
        std::vector<Face> faces;
        std::string neighbourPatch;
        std::array<double, 3> rotationAxis = {0, 0, 0};
        std::array<double, 3> rotationCenter = {0, 0, 0};
        double rotationAngle = 0;
    };
    std::map<std::string, BoundaryEntry> boundaryFacesByName;
    std::vector<Face> unnamedBoundaryFaces;
    
    for (auto& [key, info] : faceMap) {
        Face face;
        face.nodes = info.nodes;
        face.owner = info.owner;
        face.neighbour = info.neighbour;
        
        if (info.neighbour >= 0) {
            internalFaces.push_back(face);
        } else {
            if (!info.bcName.empty()) {
                auto& entry = boundaryFacesByName[info.bcName];
                entry.bcType = info.bcType;
                entry.faces.push_back(face);
                if (info.isInterface) {
                    entry.neighbourPatch = info.neighbourPatch;
                    entry.rotationAxis = info.rotationAxis;
                    entry.rotationCenter = info.rotationCenter;
                    entry.rotationAngle = info.rotationAngle;
                }
            } else {
                unnamedBoundaryFaces.push_back(face);
            }
        }
    }
    
    // Sort internal faces
    std::sort(internalFaces.begin(), internalFaces.end(), 
        [](const Face& a, const Face& b) {
            if (a.owner != b.owner) return a.owner < b.owner;
            return a.neighbour < b.neighbour;
        });
    
    // Build final face list and track face indices for faceZones
    mesh.faces.clear();
    mesh.faces.reserve(faceMap.size());
    
    // Map from faceKey to final face index
    std::unordered_map<std::string, int> faceKeyToIndex;
    
    // Add internal faces
    for (auto& f : internalFaces) {
        int faceIdx = static_cast<int>(mesh.faces.size());
        std::string key = faceKey(f.nodes);
        faceKeyToIndex[key] = faceIdx;
        mesh.faces.push_back(f);
    }
    mesh.nInternalFaces = static_cast<int>(internalFaces.size());
    
    // Add boundary faces
    mesh.boundaries.clear();
    for (auto& [bcName, entry] : boundaryFacesByName) {
        BoundaryPatch patch;
        patch.name = bcName;
        patch.type = entry.bcType;
        patch.startFace = static_cast<int>(mesh.faces.size());
        patch.nFaces = static_cast<int>(entry.faces.size());
        
        // Set cyclicAMI pairing info
        if (entry.bcType == "cyclicAMI") {
            patch.neighbourPatch = entry.neighbourPatch;
            patch.rotationAxis = entry.rotationAxis;
            patch.rotationCenter = entry.rotationCenter;
            patch.rotationAngle = entry.rotationAngle;
        }
        
        std::sort(entry.faces.begin(), entry.faces.end(),
            [](const Face& a, const Face& b) { return a.owner < b.owner; });
        
        for (auto& f : entry.faces) {
            int faceIdx = static_cast<int>(mesh.faces.size());
            std::string key = faceKey(f.nodes);
            faceKeyToIndex[key] = faceIdx;
            mesh.faces.push_back(f);
        }
        
        mesh.boundaries.push_back(patch);
    }
    
    if (!unnamedBoundaryFaces.empty()) {
        BoundaryPatch patch;
        patch.name = "defaultFaces";
        patch.type = "patch";
        patch.startFace = static_cast<int>(mesh.faces.size());
        patch.nFaces = static_cast<int>(unnamedBoundaryFaces.size());
        
        std::sort(unnamedBoundaryFaces.begin(), unnamedBoundaryFaces.end(),
            [](const Face& a, const Face& b) { return a.owner < b.owner; });
        
        for (auto& f : unnamedBoundaryFaces) {
            int faceIdx = static_cast<int>(mesh.faces.size());
            std::string key = faceKey(f.nodes);
            faceKeyToIndex[key] = faceIdx;
            mesh.faces.push_back(f);
        }
        
        mesh.boundaries.push_back(patch);
    }
    
    // Create cell zones from zone names
    for (size_t zi = 0; zi < zones.size(); ++zi) {
        Zone cellZone;
        cellZone.name = zones[zi].name;
        int nCells = static_cast<int>(zones[zi].cellNodes.size());
        for (int ci = 0; ci < nCells; ++ci) {
            cellZone.indices.push_back(ci + zones[zi].cellOffset);
        }
        mesh.cellZones.push_back(cellZone);
    }
    
    // Create face zones from boundary patches
    mesh.faceZones.clear();
    for (const auto& patch : mesh.boundaries) {
        Zone faceZone;
        faceZone.name = patch.name;
        for (int i = 0; i < patch.nFaces; ++i) {
            faceZone.indices.push_back(patch.startFace + i);
        }
        mesh.faceZones.push_back(faceZone);
    }
    
    // Remove unused points and renumber
    {
        std::set<int> usedPoints;
        for (const auto& face : mesh.faces) {
            for (int n : face.nodes) {
                usedPoints.insert(n);
            }
        }
        
        if (usedPoints.size() < mesh.points.size()) {
            std::cout << "\n  Removing " << (mesh.points.size() - usedPoints.size()) 
                      << " unused points...\n";
            
            // Create old->new point mapping
            std::vector<int> oldToNew(mesh.points.size(), -1);
            std::vector<Point3D> newPoints;
            newPoints.reserve(usedPoints.size());
            
            int newIdx = 0;
            for (int oldIdx : usedPoints) {
                oldToNew[oldIdx] = newIdx++;
                newPoints.push_back(mesh.points[oldIdx]);
            }
            
            // Update face nodes
            for (auto& face : mesh.faces) {
                for (int& n : face.nodes) {
                    n = oldToNew[n];
                }
            }
            
            mesh.points = std::move(newPoints);
        }
    }
    
    std::cout << "\nMulti-zone mesh built:\n";
    std::cout << "  Total points: " << mesh.points.size() << "\n";
    std::cout << "  Total cells: " << mesh.nCells << "\n";
    std::cout << "  Internal faces: " << mesh.nInternalFaces << "\n";
    std::cout << "  Boundary patches: " << mesh.boundaries.size() << "\n";
    for (const auto& patch : mesh.boundaries) {
        std::cout << "    - " << patch.name << " (" << patch.type << "): " << patch.nFaces << " faces\n";
    }
}

}  // namespace cgns2openfoam
