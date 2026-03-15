/**
 * Write minimal CGNS files (single zone, tetra or hex) that cgns2openfoam can read.
 * Build: same as main project, link with CGNS and HDF5.
 * Usage: write_cgns_minimal [tetra|hex] <output.cgns>
 *   Default: tetra, or output path only.
 */

#include <cgnslib.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

static int write_tetra(const char* outpath) {
    int fn, B, Z, C, S;
    if (cg_open(outpath, CG_MODE_WRITE, &fn) != CG_OK) {
        std::cerr << "cg_open write failed: " << cg_get_error() << std::endl;
        return 1;
    }
    if (cg_base_write(fn, "Base", 3, 3, &B) != CG_OK) {
        std::cerr << "cg_base_write failed: " << cg_get_error() << std::endl;
        cg_close(fn);
        return 1;
    }
    const cgsize_t nVertex = 8, nCell = 5;
    cgsize_t size[9] = { nVertex, nCell, 0, 0, 0, 0, 0, 0, 0 };
    if (cg_zone_write(fn, B, "Zone1", size, Unstructured, &Z) != CG_OK) {
        std::cerr << "cg_zone_write failed: " << cg_get_error() << std::endl;
        cg_close(fn);
        return 1;
    }
    double x[8] = {0,1,1,0, 0,1,1,0};
    double y[8] = {0,0,1,1, 0,0,1,1};
    double z[8] = {0,0,0,0, 1,1,1,1};
    if (cg_coord_write(fn, B, Z, RealDouble, "CoordinateX", x, &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateY", y, &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateZ", z, &C) != CG_OK) {
        std::cerr << "cg_coord_write failed: " << cg_get_error() << std::endl;
        cg_close(fn);
        return 1;
    }
    const cgsize_t tet_conn[5*4] = {
        1,2,3,6, 1,3,4,6, 1,4,5,6, 1,5,8,6, 1,8,2,6
    };
    if (cg_section_write(fn, B, Z, "Tetra", TETRA_4, 1, 5, 0, tet_conn, &S) != CG_OK) {
        std::cerr << "cg_section_write failed: " << cg_get_error() << std::endl;
        cg_close(fn);
        return 1;
    }
    cg_close(fn);
    std::cout << "Wrote " << outpath << " (8 points, 5 tetras)\n";
    return 0;
}

static int write_hex(const char* outpath) {
    int fn, B, Z, C, S;
    if (cg_open(outpath, CG_MODE_WRITE, &fn) != CG_OK) {
        std::cerr << "cg_open write failed: " << cg_get_error() << std::endl;
        return 1;
    }
    if (cg_base_write(fn, "Base", 3, 3, &B) != CG_OK) {
        std::cerr << "cg_base_write failed: " << cg_get_error() << std::endl;
        cg_close(fn);
        return 1;
    }
    /* 2x2x2 hex = 27 points, 8 hexes */
    const cgsize_t nVertex = 27;
    const cgsize_t nCell = 8;
    cgsize_t size[9] = { nVertex, nCell, 0, 0, 0, 0, 0, 0, 0 };
    if (cg_zone_write(fn, B, "Zone1", size, Unstructured, &Z) != CG_OK) {
        std::cerr << "cg_zone_write failed: " << cg_get_error() << std::endl;
        cg_close(fn);
        return 1;
    }
    double x[27], y[27], z[27];
    for (int iz = 0; iz < 3; iz++)
        for (int iy = 0; iy < 3; iy++)
            for (int ix = 0; ix < 3; ix++) {
                int i = ix + iy*3 + iz*9;
                x[i] = ix; y[i] = iy; z[i] = iz;
            }
    if (cg_coord_write(fn, B, Z, RealDouble, "CoordinateX", x, &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateY", y, &C) != CG_OK ||
        cg_coord_write(fn, B, Z, RealDouble, "CoordinateZ", z, &C) != CG_OK) {
        std::cerr << "cg_coord_write failed: " << cg_get_error() << std::endl;
        cg_close(fn);
        return 1;
    }
    /* CGNS hex order: 0,1,2,3 bottom, 4,5,6,7 top (1-based) */
    auto node = [](int ix, int iy, int iz) { return 1 + ix + iy*3 + iz*9; };
    cgsize_t hex_conn[8*8];
    int h = 0;
    for (int iz = 0; iz < 2; iz++)
        for (int iy = 0; iy < 2; iy++)
            for (int ix = 0; ix < 2; ix++) {
                hex_conn[h++] = node(ix,iy,iz);
                hex_conn[h++] = node(ix+1,iy,iz);
                hex_conn[h++] = node(ix+1,iy+1,iz);
                hex_conn[h++] = node(ix,iy+1,iz);
                hex_conn[h++] = node(ix,iy,iz+1);
                hex_conn[h++] = node(ix+1,iy,iz+1);
                hex_conn[h++] = node(ix+1,iy+1,iz+1);
                hex_conn[h++] = node(ix,iy+1,iz+1);
            }
    if (cg_section_write(fn, B, Z, "Hex", HEXA_8, 1, 8, 0, hex_conn, &S) != CG_OK) {
        std::cerr << "cg_section_write failed: " << cg_get_error() << std::endl;
        cg_close(fn);
        return 1;
    }
    cg_close(fn);
    std::cout << "Wrote " << outpath << " (27 points, 8 hexes)\n";
    return 0;
}

int main(int argc, char** argv) {
    const char* outpath = "minimal.cgns";
    const char* kind = "tetra";
    if (argc >= 2) {
        if (strcmp(argv[1], "tetra") == 0 || strcmp(argv[1], "hex") == 0) {
            kind = argv[1];
            outpath = argc >= 3 ? argv[2] : "minimal.cgns";
        } else {
            outpath = argv[1];
        }
    }
    if (strcmp(kind, "hex") == 0)
        return write_hex(outpath);
    return write_tetra(outpath);
}
