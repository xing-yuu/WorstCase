"""Generate an original structured beam fixture; no external model assets."""
from pathlib import Path
import argparse
import itertools
from collections import Counter


def generate(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    nx, ny, nz = 4, 3, 3
    def index(x, y, z):
        return x + (nx + 1) * (y + (ny + 1) * z)
    vertices = [(x / nx, 0.3 * y / ny, 0.2 * z / nz)
                for z in range(nz + 1) for y in range(ny + 1) for x in range(nx + 1)]
    tets = []
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                for order in itertools.permutations(range(3)):
                    v = [x, y, z]
                    tet = [index(*v)]
                    for axis in order:
                        v[axis] += 1
                        tet.append(index(*v))
                    tets.append(tet)
    with (directory / "beam.mesh").open("w", encoding="ascii", newline="\n") as out:
        out.write(f"MeshVersionFormatted 2\nDimension 3\nVertices\n{len(vertices)}\n")
        for vertex in vertices:
            out.write(" ".join(format(v, ".17g") for v in vertex) + " 0\n")
        out.write(f"Tetrahedra\n{len(tets)}\n")
        for tet in tets:
            out.write(" ".join(str(v + 1) for v in tet) + " 0\n")
        out.write("End\n")
    boundary = Counter()
    oriented = {}
    for tet in tets:
        for opposite in range(4):
            face = [tet[k] for k in range(4) if k != opposite]
            a, b, c = (vertices[v] for v in face)
            u = [b[k] - a[k] for k in range(3)]
            v = [c[k] - a[k] for k in range(3)]
            normal = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
            if sum(normal[k]*(vertices[tet[opposite]][k]-a[k]) for k in range(3)) > 0:
                face[1], face[2] = face[2], face[1]
            key = tuple(sorted(face))
            boundary[key] += 1
            oriented[key] = face
    faces = [oriented[key] for key in sorted(boundary) if boundary[key] == 1]
    used = sorted(set(v for face in faces for v in face))
    mapping = {v: i+1 for i, v in enumerate(used)}
    with (directory / "beam.obj").open("w", encoding="ascii", newline="\n") as out:
        for v in used:
            out.write("v " + " ".join(format(x, ".17g") for x in vertices[v]) + "\n")
        for face in faces:
            out.write("f " + " ".join(str(mapping[v]) for v in face) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=str(Path(__file__).resolve().parents[1] / "data"))
    generate(parser.parse_args().output)
