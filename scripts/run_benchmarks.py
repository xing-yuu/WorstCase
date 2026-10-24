"""Run supplied models and collect comparable corrected-code statistics."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True, type=Path)
    parser.add_argument("--models", required=True, nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--config", type=Path,
                        default=Path(__file__).resolve().parents[1] / "configs/paper.json")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--solver", choices=["cg", "gcr"], default="gcr")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = []
    failed = False
    for i, model in enumerate(args.models):
        directory = args.output / f"{i:03d}_{model.stem}"
        command = [str(args.exe.resolve()), "--input", str(model.resolve()),
                   "--output", str(directory.resolve()), "--config", str(args.config.resolve()),
                   "--seed", str(args.seed), "--threads", str(args.threads),
                   "--solver", args.solver]
        started = time.perf_counter()
        result = subprocess.run(command, capture_output=True, text=True)
        elapsed = time.perf_counter() - started
        (args.output / f"{directory.name}.log").write_text(
            result.stdout + "\n" + result.stderr, encoding="utf-8")
        record = {"model": model.name, "run": directory.name, "exit_code": result.returncode,
                  "wall_seconds": elapsed}
        if result.returncode == 0:
            record.update(json.loads((directory / "summary.json").read_text()))
            with (directory / "history.csv").open(newline="") as f:
                attempts = [r for r in csv.DictReader(f) if r["event"].startswith("step_")]
            record["attempted_steps"] = len(attempts)
            record["accepted_steps"] = sum(r["accepted"] == "1" for r in attempts)
        else:
            failed = True
        records.append(record)
        print(f"{model.name}: exit={result.returncode}, wall={elapsed:.3f}s", flush=True)
    (args.output / "statistics.json").write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    fields = list(dict.fromkeys(key for record in records for key in record))
    with (args.output / "statistics.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
