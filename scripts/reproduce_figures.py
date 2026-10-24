"""Plot supplied run histories; this does not reconstruct unavailable paper experiments."""
import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", required=True, type=Path)
    args = parser.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    statistics = json.loads((args.runs / "statistics.json").read_text())
    for record in statistics:
        if record["exit_code"]:
            continue
        directory = args.runs / record["run"]
        with (directory / "history.csv").open(newline="") as f:
            rows = [r for r in csv.DictReader(f)
                    if r["event"].startswith("step_") or r["event"] in ("initial_analysis", "final_analysis")]
        fig, axes = plt.subplots(1, 2, figsize=(10, 3.5), constrained_layout=True)
        x = range(len(rows))
        for key, label in [("max_boundary_stress", "Boundary tets"), ("max_all_stress", "All tets")]:
            axes[0].plot(x, [float(r[key]) for r in rows], label=label)
        axes[0].set(ylabel="Maximum absolute principal stress", xlabel="Recorded analysis / step")
        axes[0].legend()
        axes[1].plot(x, [float(r["sigma_bar"]) for r in rows], label="Barrier bound")
        axes[1].plot(x, [float(r["sigma_star"]) for r in rows], label="Target bound", linestyle="--")
        axes[1].set(xlabel="Recorded analysis / step", ylabel="Stress bound")
        axes[1].legend()
        fig.suptitle(record["model"] + " (last point recomputes the critical load)")
        for suffix in ("png", "pdf"):
            fig.savefig(directory / f"history.{suffix}", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    main()
