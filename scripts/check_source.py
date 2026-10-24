"""Check encoding and prevent private build paths from returning to release code."""
from pathlib import Path
import re
import argparse
import shutil
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--format", action="store_true", help="Also check clang-format")
options = parser.parse_args()

root = Path(__file__).resolve().parents[1]
errors = []
sources = []
for directory in ["include", "src", "apps", "tests", "scripts", "cmake"]:
    for path in (root / directory).rglob("*"):
        if path.suffix not in {".h", ".cpp", ".py", ".cmake"}:
            continue
        data = path.read_bytes()
        if path.suffix in {".h", ".cpp"}:
            sources.append(str(path))
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            errors.append(f"{path.relative_to(root)}: not UTF-8")
            continue
        if data.startswith(b"\xef\xbb\xbf") or any(ord(c) > 127 for c in text):
            errors.append(f"{path.relative_to(root)}: BOM or non-ASCII source")
        if re.search(r"(?<![A-Za-z0-9_])[A-Za-z]:[\\/]", text):
            errors.append(f"{path.relative_to(root)}: machine-specific absolute path")
        if path.suffix == ".h" and "#pragma once" not in text:
            errors.append(f"{path.relative_to(root)}: missing include guard")
        if path.suffix == ".h" and "using namespace " in text:
            errors.append(f"{path.relative_to(root)}: using namespace in public header")
if options.format:
    formatter = shutil.which("clang-format")
    if not formatter:
        errors.append("clang-format is required for --format")
    elif subprocess.run([formatter, "--dry-run", "--Werror", *sources]).returncode:
        errors.append("clang-format check failed")
if errors:
    raise SystemExit("\n".join(errors))
print("Source hygiene checks passed")
