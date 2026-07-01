"""Plot z times the collinear first moment of the pion Collins function.

For the JAM3D Gaussian Collins TMD,

    H1perp(z, pT) = 2 z^2 Mh^2 / width * collinear(z)
                    * exp(-z^2 pT^2 / width) / (pi width)

the requested moment

    H1perp(1)(z) = z^2 int d^2pT pT^2 / (2 Mh^2) H1perp(z, z^2 pT^2)

is exactly the collinear Collins factor returned by
``conf["collinspi"].get_C(z, Q2)``. This script plots z H1perp(1)(z).
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from generate_collins_grid import FLAVOR_INDEX, import_tmd, linspace, load_tmd
from generate_collins_grid import select_replicas, validate_flavors


def compute_collinear_collins(
    tmd,
    z_grid: list[float],
    q2: float,
    flavors: list[str],
    replicas: list[int],
):
    from tools.config import conf

    collins = conf["collinspi"]
    flavor_indices = [(flavor, FLAVOR_INDEX[flavor]) for flavor in flavors]
    sums = [[0.0 for _ in flavor_indices] for _ in z_grid]
    sums_sq = [[0.0 for _ in flavor_indices] for _ in z_grid]

    for irep in replicas:
        tmd.parman.set_new_params(tmd.par[irep], initial=True)
        for z_index, z in enumerate(z_grid):
            values = collins.get_C(z, q2)
            for flavor_position, (_flavor, flavor_index) in enumerate(flavor_indices):
                value = z * float(values[flavor_index])
                sums[z_index][flavor_position] += value
                sums_sq[z_index][flavor_position] += value * value

    nrep = len(replicas)
    for z_index, z in enumerate(z_grid):
        for flavor_position, (flavor, flavor_index) in enumerate(flavor_indices):
            mean = sums[z_index][flavor_position] / nrep
            mean_sq = sums_sq[z_index][flavor_position] / nrep
            variance = max(mean_sq - mean * mean, 0.0)
            yield {
                "z": z,
                "Q2": q2,
                "flavor": flavor,
                "flavor_index": flavor_index,
                "nrep": nrep,
                "zH1perp1_mean": mean,
                "zH1perp1_std": math.sqrt(variance),
            }


def write_rows(rows: list[dict[str, float | int | str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "z",
        "Q2",
        "flavor",
        "flavor_index",
        "nrep",
        "zH1perp1_mean",
        "zH1perp1_std",
    ]
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_rows(
    rows: list[dict[str, float | int | str]],
    z_grid: list[float],
    flavors: list[str],
    output: Path,
) -> None:
    ncols = 2
    nrows = math.ceil(len(flavors) / ncols)
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(9.0, 2.1 * nrows),
        sharex=True,
    )
    flat_axes = list(axes.flat)

    for ax, flavor in zip(flat_axes, flavors):
        flavor_rows = [row for row in rows if row["flavor"] == flavor]
        means = [float(row["zH1perp1_mean"]) for row in flavor_rows]
        stds = [float(row["zH1perp1_std"]) for row in flavor_rows]
        lower = [mean - std for mean, std in zip(means, stds)]
        upper = [mean + std for mean, std in zip(means, stds)]
        line = ax.plot(z_grid, means, linewidth=1.8)[0]
        ax.fill_between(z_grid, lower, upper, color=line.get_color(), alpha=0.14)
        ax.axhline(0.0, color="black", linewidth=0.8)
        ax.set_title(flavor)
        ax.grid(True, alpha=0.25)

    for ax in flat_axes[len(flavors):]:
        ax.set_visible(False)

    for ax in flat_axes[-ncols:]:
        ax.set_xlabel("z")
    for ax in flat_axes[::ncols]:
        ax.set_ylabel(r"$z H_1^{\perp(1)}(z)$")

    fig.suptitle("JAM3D pion Collins z-weighted first moment by flavor", y=0.995)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200)
    plt.close(fig)


def main() -> None:
    Q2 = 4.0
    Z_GRID = linspace(0.01, 0.99, 200)
    FLAVORS = ["u", "ub", "d", "db", "s", "sb", "c", "cb", "b", "bb"]
    REPLICA_SELECTION = 50
    TAG = "JAM3D_2022"
    CSV_OUTPUT = Path("z_collinear_collins.csv")
    PLOT_OUTPUT = Path("z_collinear_collins.png")
    JAM3DLIB_PATH = None
    JAM3D_DEV_LIB_PATH = None

    flavors = validate_flavors(FLAVORS)
    TMD, jam3dlib_path, jam3d_dev_lib_path = import_tmd(
        JAM3DLIB_PATH,
        JAM3D_DEV_LIB_PATH,
    )
    tmd = load_tmd(TMD, TAG, jam3dlib_path)
    replicas = select_replicas(REPLICA_SELECTION, tmd.nrep)
    if not replicas:
        raise SystemExit("No valid replicas selected.")

    rows = list(
        compute_collinear_collins(
            tmd=tmd,
            z_grid=Z_GRID,
            q2=Q2,
            flavors=flavors,
            replicas=replicas,
        )
    )
    write_rows(rows, CSV_OUTPUT)
    plot_rows(rows, Z_GRID, flavors, PLOT_OUTPUT)

    print(
        f"Wrote {CSV_OUTPUT} and {PLOT_OUTPUT} for {len(flavors)} flavors, "
        f"{len(Z_GRID)} z-points, and {len(replicas)} replicas.\n"
        f"jam3dlib: {jam3dlib_path}\n"
        f"jam3d_dev_lib: {jam3d_dev_lib_path}"
    )


if __name__ == "__main__":
    main()
