import gmsh

gmsh.initialize()
gmsh.model.add("cylinder")

# ---------------------
# Parameters (µm)
# ---------------------
L  = 15
D = 1
R  = D/2.0
lc = 0.2

# ---------------------
# Geometry (OpenCASCADE)
# ---------------------
gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, L, R, tag=1)
gmsh.model.occ.synchronize()

# ---------------------
# Mesh options
# ---------------------
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc)
gmsh.option.setNumber("Mesh.ElementOrder", 1)   # linear tets
gmsh.option.setNumber("Mesh.Algorithm3D", 1)

# Force Gmsh 2.2 ASCII format
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.Binary", 0)

# ---------------------
# Physical group
# ---------------------
gmsh.model.addPhysicalGroup(3, [1], name="Cylinder")

# ---------------------
# Mesh + export
# ---------------------
gmsh.model.mesh.generate(3)

# Count tetrahedra (element type 4 = 4-node tet)
elemTypes, elemTags, _ = gmsh.model.mesh.getElements(dim=3)
nTets = sum(len(tags) for et, tags in zip(elemTypes, elemTags) if et == 4)

filename = f"Cylinder_dia{D}um_L{L}um_size{lc}_{nTets}tets.msh"
gmsh.write(filename)

print(f"Exported {filename} (MSH 2.2 ASCII)")
print(f"Number of tetrahedra: {nTets}")

gmsh.finalize()
