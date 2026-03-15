#!/usr/bin/env python3
"""
Generate CGNS test meshes for cgns2openfoam test library.
Output: one mesh.cgns per case in ../01_single_zone_tetra, ../02_single_zone_hex, etc.

Requires: numpy, meshio (pip install -r requirements.txt)
CGNS written by meshio should use standard names (CoordinateX/Y/Z); if the converter
fails, ensure your meshio version writes CGNS with these coordinate names.
"""

from pathlib import Path
import sys

try:
    import numpy as np
    import meshio
except ImportError as e:
    print("Please install: pip install numpy meshio", file=sys.stderr)
    raise SystemExit(1) from e

# Case directory is sibling of scripts/
SCRIPT_DIR = Path(__file__).resolve().parent
TEST_LIB_ROOT = SCRIPT_DIR.parent
CASES_DIR = TEST_LIB_ROOT


def _ensure_cgns_coordinate_names(cgns_path: Path) -> None:
    """
    CGNS reader expects 'CoordinateX', 'CoordinateY', 'CoordinateZ'.
    meshio may write 'X'/'Y'/'Z' or 'CoordinateX' etc. Try renaming with h5py if needed.
    """
    try:
        import h5py
    except ImportError:
        return
    with h5py.File(cgns_path, "r+") as f:
        def visit(node, path=""):
            if getattr(node, "keys", None) is None:
                return
            for key in list(node.keys()):
                child = node[key]
                full = f"{path}/{key}" if path else key
                if key in ("X", "Y", "Z") and hasattr(child, "shape"):
                    new_name = f"Coordinate{key}"
                    if new_name not in node:
                        node[new_name] = child[...]
                        del node[key]
                visit(child, full)
        visit(f)


def box_tetra(Lx=1.0, Ly=1.0, Lz=1.0, nx=2, ny=2, nz=2):
    """Simple box mesh with tetrahedra (subdivide hexes into 5 tets)."""
    x = np.linspace(0, Lx, nx + 1)
    y = np.linspace(0, Ly, ny + 1)
    z = np.linspace(0, Lz, nz + 1)
    pts = []
    for zi in z:
        for yi in y:
            for xi in x:
                pts.append([xi, yi, zi])
    points = np.array(pts)
    npx, npy, npz = nx + 1, ny + 1, nz + 1

    def node_id(ix, iy, iz):
        return ix + iy * npx + iz * (npx * npy)

    tets = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                # One hex -> 5 tets (standard subdivision)
                n0 = node_id(ix, iy, iz)
                n1 = node_id(ix + 1, iy, iz)
                n2 = node_id(ix + 1, iy + 1, iz)
                n3 = node_id(ix, iy + 1, iz)
                n4 = node_id(ix, iy, iz + 1)
                n5 = node_id(ix + 1, iy, iz + 1)
                n6 = node_id(ix + 1, iy + 1, iz + 1)
                n7 = node_id(ix, iy + 1, iz + 1)
                tets.append([n0, n1, n2, n6])
                tets.append([n0, n2, n3, n6])
                tets.append([n0, n3, n7, n6])
                tets.append([n0, n7, n4, n6])
                tets.append([n0, n4, n5, n6])
    cells = [meshio.CellBlock("tetra", np.array(tets))]
    return meshio.Mesh(points, cells)


def box_hex(Lx=1.0, Ly=1.0, Lz=1.0, nx=2, ny=2, nz=2):
    """Structured hex mesh (one zone)."""
    x = np.linspace(0, Lx, nx + 1)
    y = np.linspace(0, Ly, ny + 1)
    z = np.linspace(0, Lz, nz + 1)
    pts = []
    for zi in z:
        for yi in y:
            for xi in x:
                pts.append([xi, yi, zi])
    points = np.array(pts)
    npx, npy, npz = nx + 1, ny + 1, nz + 1

    def node_id(ix, iy, iz):
        return ix + iy * npx + iz * (npx * npy)

    hexes = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                # CGNS / OpenFOAM hex order: 0,1,2,3 bottom, 4,5,6,7 top (counter-clockwise)
                n0 = node_id(ix, iy, iz)
                n1 = node_id(ix + 1, iy, iz)
                n2 = node_id(ix + 1, iy + 1, iz)
                n3 = node_id(ix, iy + 1, iz)
                n4 = node_id(ix, iy, iz + 1)
                n5 = node_id(ix + 1, iy, iz + 1)
                n6 = node_id(ix + 1, iy + 1, iz + 1)
                n7 = node_id(ix, iy + 1, iz + 1)
                hexes.append([n0, n1, n2, n3, n4, n5, n6, n7])
    cells = [meshio.CellBlock("hexahedron", np.array(hexes))]
    return meshio.Mesh(points, cells)


def channel_hex(Lx=5.0, Ly=1.0, Lz=1.0, nx=5, ny=2, nz=2):
    """Rectangular channel (duct) - hex."""
    return box_hex(Lx, Ly, Lz, nx, ny, nz)


def external_flow_box(L=2.0, n=4):
    """Cube for external flow (farfield/symmetry)."""
    return box_tetra(L, L, L, n, n, n)


def pipe_like_hex(L=3.0, D=1.0, n_axial=4, n_radial=2):
    """Approximate pipe: use a box (square cross-section) for simplicity."""
    return box_hex(L, D, D, n_axial, n_radial, n_radial)


def two_zone_interface():
    """Two hex zones sharing an interface (stator/rotor style, no rotation)."""
    # Zone 1: 0..1 in x; Zone 2: 1..2 in x. Shared face at x=1.
    pts1, hex1 = [], []
    nx, ny, nz = 2, 2, 2
    for iz in range(nz + 1):
        for iy in range(ny + 1):
            for ix in range(nx + 1):
                pts1.append([float(ix), float(iy), float(iz)])
    npts = len(pts1)
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                n0 = ix + iy * (nx + 1) + iz * (nx + 1) * (ny + 1)
                n1 = n0 + 1
                n2 = n0 + (nx + 1) + 1
                n3 = n0 + (nx + 1)
                n4 = n0 + (nx + 1) * (ny + 1)
                n5 = n1 + (nx + 1) * (ny + 1)
                n6 = n2 + (nx + 1) * (ny + 1)
                n7 = n3 + (nx + 1) * (ny + 1)
                hex1.append([n0, n1, n2, n3, n4, n5, n6, n7])
    pts2 = []
    for iz in range(nz + 1):
        for iy in range(ny + 1):
            for ix in range(nx + 1):
                pts2.append([1.0 + float(ix), float(iy), float(iz)])
    hex2 = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                n0 = ix + iy * (nx + 1) + iz * (nx + 1) * (ny + 1)
                n1 = n0 + 1
                n2 = n0 + (nx + 1) + 1
                n3 = n0 + (nx + 1)
                n4 = n0 + (nx + 1) * (ny + 1)
                n5 = n1 + (nx + 1) * (ny + 1)
                n6 = n2 + (nx + 1) * (ny + 1)
                n7 = n3 + (nx + 1) * (ny + 1)
                hex2.append([n0 + npts, n1 + npts, n2 + npts, n3 + npts,
                             n4 + npts, n5 + npts, n6 + npts, n7 + npts])
    # Single mesh with two zones: meshio doesn't support multi-zone CGNS in one Mesh easily.
    # So we output one zone (combined) for "two zone" as a single zone for now;
    # true two-zone with interface would require C++ or manual HDF5/CGNS write.
    all_pts = np.vstack([np.array(pts1), np.array(pts2)])
    all_hex = np.vstack([np.array(hex1), np.array(hex2)])
    return meshio.Mesh(all_pts, [meshio.CellBlock("hexahedron", all_hex)])


def heat_sink_like():
    """Simple block for heat sink (hex)."""
    return box_hex(1.0, 0.5, 0.3, 3, 2, 2)


def electronics_solid():
    """Solid block for temperature field."""
    return box_hex(0.1, 0.1, 0.05, 2, 2, 1)


def multi_zone_mixed():
    """Mixed: tetra + hex regions. Meshio single mesh with both cell types."""
    tet_mesh = box_tetra(0.5, 0.5, 0.5, 1, 1, 1)
    hex_mesh = box_hex(0.5, 0.5, 0.5, 1, 1, 1)
    # Offset hex to sit next to tet
    hex_mesh.points[:, 0] += 0.5
    pts = np.vstack([tet_mesh.points, hex_mesh.points])
    n_tet_pts = len(tet_mesh.points)
    hex_cells = hex_mesh.cells[0].data + n_tet_pts
    cells = [
        meshio.CellBlock("tetra", tet_mesh.cells[0].data),
        meshio.CellBlock("hexahedron", hex_cells),
    ]
    return meshio.Mesh(pts, cells)


CASE_GENERATORS = {
    "01_single_zone_tetra": lambda: box_tetra(1.0, 1.0, 1.0, 2, 2, 2),
    "02_single_zone_hex": lambda: box_hex(1.0, 1.0, 1.0, 2, 2, 2),
    "03_external_flow_3d": lambda: external_flow_box(2.0, 3),
    "04_pipe_flow": lambda: pipe_like_hex(3.0, 1.0, 4, 2),
    "05_duct_flow": lambda: channel_hex(5.0, 1.0, 1.0, 5, 2, 2),
    "06_electronics_solid": lambda: electronics_solid(),
    "07_electronics_fluid_solid": lambda: box_hex(0.2, 0.2, 0.1, 2, 2, 1),  # single zone placeholder
    "08_heat_sink": lambda: heat_sink_like(),
    "09_two_zone_interface": lambda: two_zone_interface(),
    "10_fan_rotor_stator": lambda: two_zone_interface(),  # same geometry, use for interface
    "11_fan_transient_interface": lambda: two_zone_interface(),
    "12_multi_zone_mixed": lambda: multi_zone_mixed(),
}


def main():
    for case_name, gen in CASE_GENERATORS.items():
        out_dir = CASES_DIR / case_name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / "mesh.cgns"
        print(f"Generating {case_name} -> {out_file}")
        mesh = gen()
        meshio.write(out_file, mesh, file_format="cgns")
        _ensure_cgns_coordinate_names(out_file)
    print("Done. Run cgns2openfoam on each mesh.cgns to verify.")


if __name__ == "__main__":
    main()
