/**
 * Generate CGNS test meshes for cgns_test_lib cases 01-12.
 * Usage: write_cgns_cases <case_number> <output.cgns>
 *   case_number: 01..12
 * Build: same as cgns2openfoam (CGNS + HDF5).
 */

#include <cgnslib.h>
#include <cstring>
#include <iostream>
#include <vector>

static int fail(const char* msg, int fn) {
    std::cerr << msg << cg_get_error() << std::endl;
    if (fn > 0) cg_close(fn);
    return 1;
}

/* --- 01: single zone tetra (8 pts, 5 tets) --- */
static int write_01(int fn, int B, const char* outpath) {
    int Z, C, S;
    const cgsize_t nVertex = 8, nCell = 5;
    cgsize_t size[9] = { nVertex, nCell, 0, 0, 0, 0, 0, 0, 0 };
    if (cg_zone_write(fn, B, "Zone1", size, Unstructured, &Z) != CG_OK)
        return fail("cg_zone_write: ", fn);
    double x[8] = {0,1,1,0, 0,1,1,0};
    double y[8] = {0,0,1,1, 0,0,1,1};
    double z[8] = {0,0,0,0, 1,1,1,1};
    if (cg_coord_write(fn, B, Z, RealDouble, "CoordinateX", x, &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateY", y, &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateZ", z, &C) != CG_OK)
        return fail("cg_coord_write: ", fn);
    const cgsize_t conn[20] = {1,2,3,6, 1,3,4,6, 1,4,5,6, 1,5,8,6, 1,8,2,6};
    if (cg_section_write(fn, B, Z, "Tetra", TETRA_4, 1, 5, 0, conn, &S) != CG_OK)
        return fail("cg_section_write: ", fn);
    std::cout << "Wrote " << outpath << " [01] 8 pts, 5 tetras\n";
    return 0;
}

/* --- 02: single zone hex 2x2x2 --- */
static int write_02(int fn, int B, const char* outpath) {
    int Z, C, S;
    const cgsize_t nVertex = 27, nCell = 8;
    cgsize_t size[9] = { nVertex, nCell, 0, 0, 0, 0, 0, 0, 0 };
    if (cg_zone_write(fn, B, "Zone1", size, Unstructured, &Z) != CG_OK)
        return fail("cg_zone_write: ", fn);
    double x[27], y[27], z[27];
    for (int iz = 0; iz < 3; iz++)
        for (int iy = 0; iy < 3; iy++)
            for (int ix = 0; ix < 3; ix++) {
                int i = ix + iy*3 + iz*9;
                x[i] = ix; y[i] = iy; z[i] = iz;
            }
    if (cg_coord_write(fn, B, Z, RealDouble, "CoordinateX", x, &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateY", y, &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateZ", z, &C) != CG_OK)
        return fail("cg_coord_write: ", fn);
    auto n = [](int ix, int iy, int iz) { return 1 + ix + iy*3 + iz*9; };
    cgsize_t conn[64];
    int h = 0;
    for (int iz = 0; iz < 2; iz++)
        for (int iy = 0; iy < 2; iy++)
            for (int ix = 0; ix < 2; ix++) {
                conn[h++] = n(ix,iy,iz); conn[h++] = n(ix+1,iy,iz);
                conn[h++] = n(ix+1,iy+1,iz); conn[h++] = n(ix,iy+1,iz);
                conn[h++] = n(ix,iy,iz+1); conn[h++] = n(ix+1,iy,iz+1);
                conn[h++] = n(ix+1,iy+1,iz+1); conn[h++] = n(ix,iy+1,iz+1);
            }
    if (cg_section_write(fn, B, Z, "Hex", HEXA_8, 1, 8, 0, conn, &S) != CG_OK)
        return fail("cg_section_write: ", fn);
    std::cout << "Wrote " << outpath << " [02] 27 pts, 8 hexes\n";
    return 0;
}

/* --- 03: external flow - larger tetra box (3x3x3 = 64 pts, 81 hex -> 405 tets) --- */
static int write_03(int fn, int B, const char* outpath) {
    int Z, C, S;
    const int nx = 3, ny = 3, nz = 3;
    const cgsize_t nVertex = (nx+1)*(ny+1)*(nz+1);
    const cgsize_t nHex = (cgsize_t)nx*ny*nz;
    const cgsize_t nCell = nHex * 5;
    cgsize_t size[9] = { nVertex, nCell, 0, 0, 0, 0, 0, 0, 0 };
    if (cg_zone_write(fn, B, "Zone1", size, Unstructured, &Z) != CG_OK)
        return fail("cg_zone_write: ", fn);
    std::vector<double> x(nVertex), y(nVertex), z(nVertex);
    for (int iz = 0; iz <= nz; iz++)
        for (int iy = 0; iy <= ny; iy++)
            for (int ix = 0; ix <= nx; ix++) {
                size_t i = ix + (size_t)iy*(nx+1) + (size_t)iz*(nx+1)*(ny+1);
                x[i] = ix; y[i] = iy; z[i] = iz;
            }
    if (cg_coord_write(fn, B, Z, RealDouble, "CoordinateX", x.data(), &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateY", y.data(), &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateZ", z.data(), &C) != CG_OK)
        return fail("cg_coord_write: ", fn);
    auto node = [nx,ny](int ix, int iy, int iz) {
        return (cgsize_t)(1 + ix + iy*(nx+1) + iz*(nx+1)*(ny+1));
    };
    std::vector<cgsize_t> tet_conn(nCell * 4);
    size_t idx = 0;
    for (int iz = 0; iz < nz; iz++)
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++) {
                cgsize_t n0 = node(ix,iy,iz), n1 = node(ix+1,iy,iz), n2 = node(ix+1,iy+1,iz), n3 = node(ix,iy+1,iz);
                cgsize_t n4 = node(ix,iy,iz+1), n5 = node(ix+1,iy,iz+1), n6 = node(ix+1,iy+1,iz+1), n7 = node(ix,iy+1,iz+1);
                tet_conn[idx++] = n0; tet_conn[idx++] = n1; tet_conn[idx++] = n2; tet_conn[idx++] = n6;
                tet_conn[idx++] = n0; tet_conn[idx++] = n2; tet_conn[idx++] = n3; tet_conn[idx++] = n6;
                tet_conn[idx++] = n0; tet_conn[idx++] = n3; tet_conn[idx++] = n7; tet_conn[idx++] = n6;
                tet_conn[idx++] = n0; tet_conn[idx++] = n7; tet_conn[idx++] = n4; tet_conn[idx++] = n6;
                tet_conn[idx++] = n0; tet_conn[idx++] = n4; tet_conn[idx++] = n5; tet_conn[idx++] = n6;
            }
    if (cg_section_write(fn, B, Z, "Tetra", TETRA_4, 1, (cgsize_t)nCell, 0, tet_conn.data(), &S) != CG_OK)
        return fail("cg_section_write: ", fn);
    std::cout << "Wrote " << outpath << " [03] " << nVertex << " pts, " << nCell << " tetras\n";
    return 0;
}

/* --- 04: pipe flow - hex 4x2x2 --- */
static int write_hex_box(int fn, int B, const char* outpath, int nx, int ny, int nz, const char* label) {
    int Z, C, S;
    const cgsize_t nVertex = (nx+1)*(ny+1)*(nz+1);
    const cgsize_t nCell = (cgsize_t)nx*ny*nz;
    cgsize_t size[9] = { nVertex, nCell, 0, 0, 0, 0, 0, 0, 0 };
    if (cg_zone_write(fn, B, "Zone1", size, Unstructured, &Z) != CG_OK)
        return fail("cg_zone_write: ", fn);
    std::vector<double> x(nVertex), y(nVertex), z(nVertex);
    for (int iz = 0; iz <= nz; iz++)
        for (int iy = 0; iy <= ny; iy++)
            for (int ix = 0; ix <= nx; ix++) {
                size_t i = ix + (size_t)iy*(nx+1) + (size_t)iz*(nx+1)*(ny+1);
                x[i] = ix; y[i] = iy; z[i] = iz;
            }
    if (cg_coord_write(fn, B, Z, RealDouble, "CoordinateX", x.data(), &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateY", y.data(), &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateZ", z.data(), &C) != CG_OK)
        return fail("cg_coord_write: ", fn);
    auto n = [nx,ny](int ix, int iy, int iz) {
        return (cgsize_t)(1 + ix + iy*(nx+1) + iz*(nx+1)*(ny+1));
    };
    std::vector<cgsize_t> conn(nCell * 8);
    size_t h = 0;
    for (int iz = 0; iz < nz; iz++)
        for (int iy = 0; iy < ny; iy++)
            for (int ix = 0; ix < nx; ix++) {
                conn[h++] = n(ix,iy,iz); conn[h++] = n(ix+1,iy,iz);
                conn[h++] = n(ix+1,iy+1,iz); conn[h++] = n(ix,iy+1,iz);
                conn[h++] = n(ix,iy,iz+1); conn[h++] = n(ix+1,iy,iz+1);
                conn[h++] = n(ix+1,iy+1,iz+1); conn[h++] = n(ix,iy+1,iz+1);
            }
    if (cg_section_write(fn, B, Z, "Hex", HEXA_8, 1, nCell, 0, conn.data(), &S) != CG_OK)
        return fail("cg_section_write: ", fn);
    std::cout << "Wrote " << outpath << " [" << label << "] " << nVertex << " pts, " << nCell << " hexes\n";
    return 0;
}

static int write_04(int fn, int B, const char* outpath) { return write_hex_box(fn, B, outpath, 4, 2, 2, "04"); }
static int write_05(int fn, int B, const char* outpath) { return write_hex_box(fn, B, outpath, 5, 2, 2, "05"); }
static int write_06(int fn, int B, const char* outpath) { return write_hex_box(fn, B, outpath, 2, 2, 1, "06"); }
static int write_07(int fn, int B, const char* outpath) { return write_hex_box(fn, B, outpath, 3, 3, 2, "07"); }
static int write_08(int fn, int B, const char* outpath) { return write_hex_box(fn, B, outpath, 3, 2, 2, "08"); }

/* --- 09/10/11: two zones (stator + rotor), no 1-to-1 for simplicity --- */
static int write_two_zones(int fn, int B, const char* outpath, const char* label) {
    int Z1, Z2, C, S;
    const int nx = 2, ny = 2, nz = 2;
    const cgsize_t nV = (nx+1)*(ny+1)*(nz+1);
    const cgsize_t nC = (cgsize_t)nx*ny*nz;
    cgsize_t size1[9] = { nV, nC, 0, 0, 0, 0, 0, 0, 0 };
    cgsize_t size2[9] = { nV, nC, 0, 0, 0, 0, 0, 0, 0 };
    if (cg_zone_write(fn, B, "Stator", size1, Unstructured, &Z1) != CG_OK ||
        cg_zone_write(fn, B, "Rotor",  size2, Unstructured, &Z2) != CG_OK)
        return fail("cg_zone_write: ", fn);
    double x[27], y[27], z[27];
    for (int iz = 0; iz < 3; iz++)
        for (int iy = 0; iy < 3; iy++)
            for (int ix = 0; ix < 3; ix++) {
                int i = ix + iy*3 + iz*9;
                x[i] = ix; y[i] = iy; z[i] = iz;
            }
    if (cg_coord_write(fn, B, Z1, RealDouble, "CoordinateX", x, &C) != CG_OK ||
        cg_coord_write(fn, B, Z1, RealDouble, "CoordinateY", y, &C) != CG_OK ||
        cg_coord_write(fn, B, Z1, RealDouble, "CoordinateZ", z, &C) != CG_OK)
        return fail("Zone1 coords: ", fn);
    for (int i = 0; i < 27; i++) x[i] += 2.0;
    if (cg_coord_write(fn, B, Z2, RealDouble, "CoordinateX", x, &C) != CG_OK ||
        cg_coord_write(fn, B, Z2, RealDouble, "CoordinateY", y, &C) != CG_OK ||
        cg_coord_write(fn, B, Z2, RealDouble, "CoordinateZ", z, &C) != CG_OK)
        return fail("Zone2 coords: ", fn);
    auto n = [](int ix, int iy, int iz) { return (cgsize_t)(1 + ix + iy*3 + iz*9); };
    cgsize_t conn[64];
    int h = 0;
    for (int iz = 0; iz < 2; iz++)
        for (int iy = 0; iy < 2; iy++)
            for (int ix = 0; ix < 2; ix++) {
                conn[h++] = n(ix,iy,iz); conn[h++] = n(ix+1,iy,iz);
                conn[h++] = n(ix+1,iy+1,iz); conn[h++] = n(ix,iy+1,iz);
                conn[h++] = n(ix,iy,iz+1); conn[h++] = n(ix+1,iy,iz+1);
                conn[h++] = n(ix+1,iy+1,iz+1); conn[h++] = n(ix,iy+1,iz+1);
            }
    if (cg_section_write(fn, B, Z1, "Hex", HEXA_8, 1, 8, 0, conn, &S) != CG_OK ||
        cg_section_write(fn, B, Z2, "Hex", HEXA_8, 1, 8, 0, conn, &S) != CG_OK)
        return fail("cg_section_write: ", fn);
    std::cout << "Wrote " << outpath << " [" << label << "] 2 zones, 27+27 pts, 8+8 hexes\n";
    return 0;
}

/* --- 12: one zone, mixed tetra + hex (two sections) --- */
static int write_12(int fn, int B, const char* outpath) {
    int Z, C, S;
    /* 8 pts for tetras + 27 for hexes, but we share: one block 0..1, one 1..2. Shared face at x=1: 4 points.
     * Simpler: 8 + 19 = 27 points (tet block 0-1, hex block 1-2 with shared face 4 pts). So 8+27-4 = 31 points?
     * Easiest: one zone with 35 points - 8 for cube (5 tets) + 27 for 2x2x2 hex offset. No shared: 8+27=35. */
    const cgsize_t nVertex = 35;
    const cgsize_t nCell = 5 + 8;
    cgsize_t size[9] = { nVertex, nCell, 0, 0, 0, 0, 0, 0, 0 };
    if (cg_zone_write(fn, B, "Zone1", size, Unstructured, &Z) != CG_OK)
        return fail("cg_zone_write: ", fn);
    std::vector<double> x(nVertex), y(nVertex), z(nVertex);
    for (int i = 0; i < 8; i++) {
        x[i] = (i%4 < 2) ? (i/4) : (i/4);
        y[i] = (i%2) ? 1 : 0;
        z[i] = (i/2)%2;
    }
    x[0]=0;x[1]=1;x[2]=1;x[3]=0; x[4]=0;x[5]=1;x[6]=1;x[7]=0;
    y[0]=0;y[1]=0;y[2]=1;y[3]=1; y[4]=0;y[5]=0;y[6]=1;y[7]=1;
    z[0]=0;z[1]=0;z[2]=0;z[3]=0; z[4]=1;z[5]=1;z[6]=1;z[7]=1;
    for (int iz = 0; iz < 3; iz++)
        for (int iy = 0; iy < 3; iy++)
            for (int ix = 0; ix < 3; ix++) {
                size_t i = 8 + ix + iy*3 + iz*9;
                x[i] = 1.0 + ix; y[i] = iy; z[i] = iz;
            }
    if (cg_coord_write(fn, B, Z, RealDouble, "CoordinateX", x.data(), &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateY", y.data(), &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateZ", z.data(), &C) != CG_OK)
        return fail("cg_coord_write: ", fn);
    const cgsize_t tet_conn[20] = {1,2,3,6, 1,3,4,6, 1,4,5,6, 1,5,8,6, 1,8,2,6};
    if (cg_section_write(fn, B, Z, "Tetra", TETRA_4, 1, 5, 0, tet_conn, &S) != CG_OK)
        return fail("cg_section_write Tetra: ", fn);
    cgsize_t hex_conn[64];
    int h = 0;
    for (int iz = 0; iz < 2; iz++)
        for (int iy = 0; iy < 2; iy++)
            for (int ix = 0; ix < 2; ix++) {
                cgsize_t o = 9;
                hex_conn[h++] = o + ix+iy*3+iz*9;
                hex_conn[h++] = o + ix+1+iy*3+iz*9;
                hex_conn[h++] = o + ix+1+(iy+1)*3+iz*9;
                hex_conn[h++] = o + ix+(iy+1)*3+iz*9;
                hex_conn[h++] = o + ix+iy*3+(iz+1)*9;
                hex_conn[h++] = o + ix+1+iy*3+(iz+1)*9;
                hex_conn[h++] = o + ix+1+(iy+1)*3+(iz+1)*9;
                hex_conn[h++] = o + ix+(iy+1)*3+(iz+1)*9;
            }
    if (cg_section_write(fn, B, Z, "Hex", HEXA_8, 6, 13, 0, hex_conn, &S) != CG_OK)
        return fail("cg_section_write Hex: ", fn);
    std::cout << "Wrote " << outpath << " [12] 35 pts, 5 tetras + 8 hexes\n";
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: write_cgns_cases <01|02|...|12> <output.cgns>\n";
        return 1;
    }
    int caseNum = 0;
    if (sscanf(argv[1], "%d", &caseNum) != 1 || caseNum < 1 || caseNum > 12) {
        std::cerr << "Case must be 01-12\n";
        return 1;
    }
    const char* outpath = argv[2];
    int fn, B;
    if (cg_open(outpath, CG_MODE_WRITE, &fn) != CG_OK)
        return fail("cg_open: ", 0);
    if (cg_base_write(fn, "Base", 3, 3, &B) != CG_OK)
        return fail("cg_base_write: ", fn);

    int ret = 0;
    switch (caseNum) {
        case 1: ret = write_01(fn, B, outpath); break;
        case 2: ret = write_02(fn, B, outpath); break;
        case 3: ret = write_03(fn, B, outpath); break;
        case 4: ret = write_04(fn, B, outpath); break;
        case 5: ret = write_05(fn, B, outpath); break;
        case 6: ret = write_06(fn, B, outpath); break;
        case 7: ret = write_07(fn, B, outpath); break;
        case 8: ret = write_08(fn, B, outpath); break;
        case 9: ret = write_two_zones(fn, B, outpath, "09"); break;
        case 10: ret = write_two_zones(fn, B, outpath, "10"); break;
        case 11: ret = write_two_zones(fn, B, outpath, "11"); break;
        case 12: ret = write_12(fn, B, outpath); break;
        default: ret = 1; break;
    }
    cg_close(fn);
    return ret;
}
