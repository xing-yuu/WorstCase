"""Exercise all advertised surface formats through the optional mesher."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

exe, source, output = map(Path, sys.argv[1:4])
output.mkdir(parents=True, exist_ok=True)
root = Path(tempfile.mkdtemp(prefix="surface-", dir=output))
vertices, faces = [], []
for line in (source / "data/beam.obj").read_text().splitlines():
    row = line.split()
    if row[0] == "v":
        vertices.append(list(map(float, row[1:])))
    elif row[0] == "f":
        faces.append([int(v)-1 for v in row[1:]])
ply = root / "beam.ply"
ply.write_text("ply\nformat ascii 1.0\nelement vertex " + str(len(vertices)) +
               "\nproperty double x\nproperty double y\nproperty double z\nelement face " +
               str(len(faces)) + "\nproperty list uchar int vertex_indices\nend_header\n" +
               "".join(" ".join(map(str, v)) + "\n" for v in vertices) +
               "".join("3 " + " ".join(map(str, f)) + "\n" for f in faces))
stl = root / "beam.stl"
with stl.open("w") as f:
    f.write("solid beam\n")
    for face in faces:
        f.write("facet normal 0 0 0\nouter loop\n")
        for v in face:
            f.write("vertex " + " ".join(map(str, vertices[v])) + "\n")
        f.write("endloop\nendfacet\n")
    f.write("endsolid beam\n")
config = root / "config.json"
config.write_text(json.dumps({"tet_edge_length": .25, "threads": 1, "seed": 17}))
for model in (source / "data/beam.obj", ply, stl):
    destination = root / model.suffix[1:]
    result = subprocess.run([str(exe), "--input", str(model), "--output", str(destination),
                             "--config", str(config), "--analyze-only"],
                            capture_output=True, text=True, timeout=300)
    (root / (model.suffix[1:] + ".log")).write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stderr
    report = json.loads((destination / "analysis.json").read_text())
    assert 0 < report["boundary_tets"] <= report["all_tets"]
    assert report["max_all_stress"] > 0
print(f"OBJ/PLY/STL tetrahedralization and analysis passed: {root}")
repeat = subprocess.run([str(exe), "--input", str(source / "data/beam.obj"),
                          "--output", str(root / "repeat"), "--config", str(config), "--analyze-only"],
                         capture_output=True, text=True, timeout=300)
assert repeat.returncode == 0, repeat.stderr
assert ((root / "repeat/generated_input.mesh").read_bytes() ==
        (root / "obj/generated_input.mesh").read_bytes()), "Serial meshing changed between fresh processes"
