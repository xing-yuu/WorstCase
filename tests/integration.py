"""CLI, deterministic resume, both solvers, normalization, and failure tests."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

exe, source, output = map(Path, sys.argv[1:4])
output.mkdir(parents=True, exist_ok=True)
root = Path(tempfile.mkdtemp(prefix="integration-", dir=output))
model = source / "data/beam.mesh"
before = hashlib.sha256(model.read_bytes()).hexdigest()


def run(args, success=True):
    result = subprocess.run([str(exe), *map(str, args)], cwd=root, capture_output=True,
                            text=True, timeout=300)
    if (result.returncode == 0) != success:
        raise AssertionError(f"Exit {result.returncode}: {args}\n{result.stdout}\n{result.stderr}")
    return result


run(["--help"])
config = {"outer_iterations": 2, "inner_iterations": 2, "sampling_density": 100,
          "threads": 1, "seed": 17, "remesh": False, "save_every": 1}
config_file = root / "config.json"
config_file.write_text(json.dumps(config))
for solver in ["cg", "gcr"]:
    run(["--input", model, "--output", root / solver, "--config", config_file,
         "--solver", solver, "--analyze-only"])
a = json.loads((root / "cg/analysis.json").read_text())
b = json.loads((root / "gcr/analysis.json").read_text())
assert len(a["eigenvalues"]) == len(b["eigenvalues"]) == 3
assert a["boundary_tets"] < a["all_tets"]
assert abs(a["max_boundary_stress"] - b["max_boundary_stress"]) / b["max_boundary_stress"] < 1e-5
assert a["penalty"] / a["lambda_max"] == 1e7
run(["--input", model, "--output", root / "parallel", "--config", config_file,
     "--threads", 2, "--analyze-only"])
parallel = json.loads((root / "parallel/analysis.json").read_text())
assert abs(parallel["max_boundary_stress"] / b["max_boundary_stress"] - 1) < 1e-7
run(["--input", model, "--output", root / "full", "--config", config_file])
run(["--input", model, "--output", root / "first", "--config", config_file, "--outer-iterations", 1])
run(["--resume", root / "first/checkpoint.json", "--output", root / "resumed", "--outer-iterations", 2])
assert (root / "full/final.obj").read_bytes() == (root / "resumed/final.obj").read_bytes(), "Resume changed final geometry"
for name in ["full", "first", "resumed"]:
    for filename in ["final.obj", "final.mesh", "history.csv", "config_used.json", "checkpoint.json", "summary.json"]:
        assert (root / name / filename).is_file()
rows = list(csv.DictReader((root / "full/history.csv").open()))
assert rows[0]["event"] == "initial_analysis"
assert abs(float(rows[0]["max_boundary_stress"]) / b["max_boundary_stress"] - 1) < 1e-10
assert len([r for r in rows if r["event"].startswith("step_")]) == 4
assert rows[-1]["event"] == "final_analysis"
assert {r["phase"] for r in rows if r["event"].startswith("step_")} == {"1", "2"}
for previous, current in zip(rows, rows[1:]):
    if (current["event"] == "step_accepted" and current["phase"] == previous["phase"] and
            current["outer"] == previous["outer"]):
        assert float(current["F"]) <= float(previous["F"]) + 1e-7
assert hashlib.sha256(model.read_bytes()).hexdigest() == before, "Input modified"
run(["--input", model, "--output", root / "full"], False)
run(["--input", model, "--output", root / "bad", "--solver", "unknown"], False)
bad_config = root / "bad_config.json"
for setting in [{"lame_mu": -1}, {"lame_lambda": -1, "lame_mu": .3}, {"unknown": 1}, {"threads": 1.5}]:
    bad_config.write_text(json.dumps(setting))
    run(["--input", model, "--output", root / "bad", "--config", bad_config], False)
run(["--input", model, "--output", root / "bad", "--threads", "2junk"], False)
malformed = json.loads((root / "first/checkpoint.json").read_text())
malformed["pressure"] = [[1.0]]
(root / "malformed.json").write_text(json.dumps(malformed))
run(["--resume", root / "malformed.json", "--output", root / "bad"], False)

# Exercise node-tag remapping in both supported Gmsh versions and the
# longest-side normalization using the same tetrahedra in different units.
tokens = model.read_text().split()
start = tokens.index("Vertices") + 1
count = int(tokens[start])
vertices = [[float(v) for v in tokens[start+1+4*i:start+4+4*i]] for i in range(count)]
start = tokens.index("Tetrahedra") + 1
tet_count = int(tokens[start])
tets = [[int(v)-1 for v in tokens[start+1+5*i:start+5+5*i]] for i in range(tet_count)]
tags = [100+3*i for i in range(count)]
node_lines = [" ".join(map(str, [7*v[0]+2, 7*v[1]-3, 7*v[2]+5])) for v in vertices]
for version in ("2.2", "4.1"):
    path = root / f"beam_{version}.msh"
    text = f"$MeshFormat\n{version} 0 8\n$EndMeshFormat\n$Nodes\n"
    if version == "2.2":
        text += str(count) + "\n" + "".join(f"{tags[i]} {line}\n" for i, line in enumerate(node_lines))
        text += "$EndNodes\n$Elements\n" + str(tet_count) + "\n"
        text += "".join(f"{i+1} 4 0 " + " ".join(str(tags[v]) for v in tet) + "\n"
                        for i, tet in enumerate(tets))
    else:
        text += f"1 {count} {tags[0]} {tags[-1]}\n3 1 0 {count}\n"
        text += "\n".join(map(str, tags)) + "\n" + "\n".join(node_lines) + "\n"
        text += f"$EndNodes\n$Elements\n1 {tet_count} 1 {tet_count}\n3 1 4 {tet_count}\n"
        text += "".join(f"{i+1} " + " ".join(str(tags[v]) for v in tet) + "\n"
                        for i, tet in enumerate(tets))
    path.write_text(text + "$EndElements\n")
    run(["--input", path, "--output", root / version, "--config", config_file, "--analyze-only"])
    analysis = json.loads((root / version / "analysis.json").read_text())
    assert abs(analysis["max_boundary_stress"] / b["max_boundary_stress"] - 1) < 1e-5

annotated = root / "annotated.mesh"
annotated.write_text(model.read_text().replace("End", "Corners\n1\n1\nRequiredVertices\n1\n2\nEnd"))
run(["--input", annotated, "--output", root / "annotated", "--config", config_file, "--analyze-only"])
reference_file = source / "tests/reference/beam.json"
if reference_file.is_file():
    reference = json.loads(reference_file.read_text())
    for key in ("max_boundary_stress", "max_all_stress", "sigma_bar"):
        assert abs(b[key]/reference["analysis"][key]-1) < 2e-4, f"Regression changed {key}"
    summary = json.loads((root / "full/summary.json").read_text())
    for key in ("max_boundary_stress", "max_all_stress", "E_hd"):
        assert abs(summary[key]/reference["summary"][key]-1) < 5e-3, f"Trajectory regression changed {key}"
    assert len(rows) == len(reference["history"])
    for current, expected in zip(rows, reference["history"]):
        for key in ("outer", "inner", "phase", "accepted", "event"):
            assert current[key] == expected[key], f"Scheduler regression changed {key}"
        for key in ("sigma_bar", "sigma_star", "eta", "F", "E_hd", "E_lds"):
            assert abs(float(current[key])-float(expected[key])) <= 1e-8 + 5e-3*abs(float(expected[key])), key
print(f"Integration tests passed; evidence: {root}")
