"""Plot z times the collinear first moment of the pion Collins function."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Collins import linspace, load_collins_fit, output_path, write_rows
from Collins import z_collinear_collins_rows


def plot_rows(
    rows: list[dict[str, float | int | str]],
    z_grid: list[float],
    flavors: list[str],
    output: Path,
) -> None:
    """Plot replica mean bands for z times the collinear Collins function."""
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
        means = [float(row["mean"]) for row in flavor_rows]
        stds = [float(row["std"]) for row in flavor_rows]
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
    """Write a CSV and plot of z-weighted collinear pion Collins moments."""
    Q2 = 4.0
    Z_GRID = linspace(0.20, 0.99, 200)
    FLAVORS = ["u", "ub", "d", "db", "s", "sb", "c", "cb", "b", "bb"]
    REPLICA_SELECTION = 50
    TAG = "JAM3D_2022"
    CSV_OUTPUT = output_path("z_collinear_collins.csv")
    PLOT_OUTPUT = output_path("z_collinear_collins.png")

    tmd, replicas, jam3dlib_path, jam3d_dev_lib_path = load_collins_fit(
        tag=TAG,
        replica_selection=REPLICA_SELECTION,
    )
    rows = z_collinear_collins_rows(
        tmd=tmd,
        z_grid=Z_GRID,
        q2=Q2,
        flavors=FLAVORS,
        replicas=replicas,
    )
    output_rows = [
        {
            "z": row["z"],
            "Q2": row["Q2"],
            "flavor": row["flavor"],
            "flavor_index": row["flavor_index"],
            "nrep": row["nrep"],
            "zH1perp1_mean": row["mean"],
            "zH1perp1_std": row["std"],
        }
        for row in rows
    ]
    write_rows(
        rows=output_rows,
        output=CSV_OUTPUT,
        fieldnames=[
            "z",
            "Q2",
            "flavor",
            "flavor_index",
            "nrep",
            "zH1perp1_mean",
            "zH1perp1_std",
        ],
    )
    plot_rows(rows, Z_GRID, FLAVORS, PLOT_OUTPUT)

    print(
        f"Wrote {CSV_OUTPUT} and {PLOT_OUTPUT} for {len(FLAVORS)} flavors, "
        f"{len(Z_GRID)} z-points, and {len(replicas)} replicas.\n"
        f"jam3dlib: {jam3dlib_path}\n"
        f"jam3d_dev_lib: {jam3d_dev_lib_path}"
    )


if __name__ == "__main__":
    main()
