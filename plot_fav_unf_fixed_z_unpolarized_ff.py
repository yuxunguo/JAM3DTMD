"""Plot favored and unfavored unpolarized pion TMD FF bands at fixed z and Q2."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Collins import fixed_z_unpolarized_ff_rows, linspace, load_collins_fit
from Collins import output_path, write_rows


def plot_rows(
    rows: list[dict[str, float | int | str]],
    pt_grid: list[float],
    z: float,
    q2: float,
    flavors: list[str],
    output: Path,
) -> None:
    """Plot favored/unfavored replica bands for unpolarized TMD FF values."""
    colors = {"fav": "red", "unf": "blue"}
    fig, ax = plt.subplots(figsize=(6.8, 5.0), facecolor="white")
    ax.set_facecolor("white")

    for flavor in flavors:
        flavor_rows = [row for row in rows if row["flavor"] == flavor]
        means = [float(row["mean"]) for row in flavor_rows]
        stds = [float(row["std"]) for row in flavor_rows]
        lower = [mean - std for mean, std in zip(means, stds)]
        upper = [mean + std for mean, std in zip(means, stds)]
        color = colors[flavor]
        ax.fill_between(pt_grid, lower, upper, color=color, alpha=0.28, label=flavor)
        ax.plot(pt_grid, means, color=color, linewidth=2.0)

    ax.text(
        0.98,
        0.78,
        rf"$z = {z:.1f},\ Q^2 = {q2:.0f}\ \mathrm{{GeV}}^2$",
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=10,
    )
    ax.legend(loc="upper right", frameon=True, fontsize=9)
    ax.set_xlim(min(pt_grid), max(pt_grid))
    ax.set_ylim(bottom=0.0)
    ax.set_title("JAM3D unpolarized pion TMD FF at fixed z", pad=10)
    ax.set_xlabel(r"$p_T\ \mathrm{(GeV)}$")
    ax.set_ylabel(r"$D_1(z,z p_T)$")
    ax.grid(True, which="major", color="0.75", linewidth=0.8, alpha=0.45)
    ax.minorticks_on()
    ax.grid(True, which="minor", color="0.85", linewidth=0.5, alpha=0.25)
    ax.tick_params(axis="both", which="major", direction="out", length=5, width=1)
    ax.tick_params(axis="both", which="minor", direction="out", length=3, width=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
        spine.set_color("black")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200, facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    """Write a CSV and plot of unpolarized TMD FF bands at fixed z."""
    Z = 0.3
    Q2 = 4.0
    PT_GRID = linspace(0.0, 1.0, 200)
    FLAVORS = ["fav", "unf"]
    REPLICA_SELECTION = 50
    TAG = "JAM3D_2022"
    CSV_OUTPUT = output_path("fav_unf_fixed_z_unpolarized_ff.csv")
    PLOT_OUTPUT = output_path("fav_unf_fixed_z_unpolarized_ff.png")

    tmd, replicas, jam3dlib_path, jam3d_dev_lib_path = load_collins_fit(
        tag=TAG,
        replica_selection=REPLICA_SELECTION,
    )
    rows = fixed_z_unpolarized_ff_rows(
        tmd=tmd,
        z=Z,
        pt_grid=PT_GRID,
        q2=Q2,
        flavors=FLAVORS,
        replicas=replicas,
    )
    output_rows = [
        {
            "z": row["z"],
            "pT": row["pT"],
            "z_pT": row["z_pT"],
            "Q2": row["Q2"],
            "flavor": row["flavor"],
            "flavor_index": row["flavor_index"],
            "nrep": row["nrep"],
            "unpolarized_ff_mean": row["mean"],
            "unpolarized_ff_std": row["std"],
        }
        for row in rows
    ]
    write_rows(
        rows=output_rows,
        output=CSV_OUTPUT,
        fieldnames=[
            "z",
            "pT",
            "z_pT",
            "Q2",
            "flavor",
            "flavor_index",
            "nrep",
            "unpolarized_ff_mean",
            "unpolarized_ff_std",
        ],
    )
    plot_rows(rows, PT_GRID, Z, Q2, FLAVORS, PLOT_OUTPUT)

    print(
        f"Wrote {CSV_OUTPUT} and {PLOT_OUTPUT} at z={Z}, Q2={Q2} "
        f"for {len(replicas)} replicas.\n"
        f"jam3dlib: {jam3dlib_path}\n"
        f"jam3d_dev_lib: {jam3d_dev_lib_path}"
    )


if __name__ == "__main__":
    main()
