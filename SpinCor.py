"""Calculate and plot the Collins spin-correlation observable with JAM3D TMD FFs."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from Collins import FLAVOR_INDEX, MPI_GEV, linspace, load_collins_fit
from Collins import output_path, write_rows


def mean_and_std(values: Sequence[float]) -> tuple[float, float]:
    """Return population mean and standard deviation for replica samples."""
    if not values:
        raise ValueError("cannot summarize an empty value list")
    mean = sum(values) / len(values)
    variance = max(sum(value * value for value in values) / len(values) - mean * mean, 0.0)
    return mean, math.sqrt(variance)


def tmd_ff_values(tmd, z: float, pperp: float, q2: float, irep: int) -> tuple:
    """Return JAM3D unpolarized and Collins pion TMD FF arrays.

    The equations use P_hperp = z k_perp, while the JAM3D wrapper argument is
    kT.  Therefore the wrapper is evaluated at kT = P_hperp / z.
    """
    if z <= 0.0:
        raise ValueError("z must be positive")
    kt = pperp / z
    d1 = tmd.eval(z, q2, kt, "pi", "ff", irep, icol=False)
    h1perp = tmd.eval(z, q2, kt, "pi", "collinspi", irep, icol=False)
    return d1, h1perp


def spin_correlation_rows(
    tmd,
    pperp_grid: list[float],
    z: float,
    zbar: float,
    q2: float,
    channels: list[str],
    replicas: list[int],
    cnn: float,
    mh: float = MPI_GEV,
) -> list[dict[str, float | int | str]]:
    """Build replica mean/std rows for <sin(phi) sin(phi')>.

    The reported coefficient is the JAM3D FF factor multiplying Cnn:

        coefficient = (1/4) P P' H Hbar / (z zbar mh^2 D Dbar)

    The spin_correlation column is coefficient * Cnn.  This script uses a
    symmetric transverse-momentum pair, P' = P.
    """
    rows: list[dict[str, float | int | str]] = []
    channel_indices = [(channel, FLAVOR_INDEX[channel]) for channel in channels]

    for pperp in pperp_grid:
        pperp_bar = pperp
        samples = {channel: [] for channel in channels}
        coefficient_samples = {channel: [] for channel in channels}
        d1_samples = {channel: [] for channel in channels}
        h1_samples = {channel: [] for channel in channels}

        for irep in replicas:
            d1, h1perp = tmd_ff_values(tmd, z, pperp, q2, irep)
            d1_bar, h1perp_bar = tmd_ff_values(tmd, zbar, pperp_bar, q2, irep)

            for channel, flavor_index in channel_indices:
                d_q = float(d1[flavor_index])
                d_qbar = float(d1_bar[flavor_index])
                h_q = float(h1perp[flavor_index])
                h_qbar = float(h1perp_bar[flavor_index])

                coefficient = (
                    0.25
                    * pperp
                    * pperp_bar
                    * h_q
                    * h_qbar
                    / (z * zbar * mh * mh * d_q * d_qbar)
                )
                coefficient_samples[channel].append(coefficient)
                samples[channel].append(cnn * coefficient)
                d1_samples[channel].append(d_q)
                h1_samples[channel].append(h_q)

        for channel, flavor_index in channel_indices:
            mean, std = mean_and_std(samples[channel])
            coefficient_mean, coefficient_std = mean_and_std(
                coefficient_samples[channel]
            )
            d1_mean, d1_std = mean_and_std(d1_samples[channel])
            h1_mean, h1_std = mean_and_std(h1_samples[channel])
            rows.append(
                {
                    "P_perp": pperp,
                    "Pbar_perp": pperp_bar,
                    "z": z,
                    "zbar": zbar,
                    "Q2": q2,
                    "Cnn": cnn,
                    "mh": mh,
                    "channel": channel,
                    "flavor_index": flavor_index,
                    "nrep": len(replicas),
                    "coefficient_mean": coefficient_mean,
                    "coefficient_std": coefficient_std,
                    "spin_correlation_mean": mean,
                    "spin_correlation_std": std,
                    "D1_mean": d1_mean,
                    "D1_std": d1_std,
                    "H1perp_mean": h1_mean,
                    "H1perp_std": h1_std,
                }
            )

    return rows


def plot_rows(
    rows: list[dict[str, float | int | str]],
    pperp_grid: list[float],
    channels: list[str],
    z: float,
    zbar: float,
    q2: float,
    cnn: float,
    output: Path,
) -> None:
    """Plot the extracted spin correlation for each favored/unfavored channel."""
    colors = {"fav": "red", "unf": "blue"}
    fig, ax = plt.subplots(figsize=(7.0, 5.0), facecolor="white")
    ax.set_facecolor("white")

    for channel in channels:
        channel_rows = [row for row in rows if row["channel"] == channel]
        means = [float(row["spin_correlation_mean"]) for row in channel_rows]
        stds = [float(row["spin_correlation_std"]) for row in channel_rows]
        lower = [mean - std for mean, std in zip(means, stds)]
        upper = [mean + std for mean, std in zip(means, stds)]
        color = colors[channel]
        ax.fill_between(pperp_grid, lower, upper, color=color, alpha=0.25)
        ax.plot(pperp_grid, means, color=color, linewidth=2.0, label=channel)

    ax.axhline(0.0, color="black", linewidth=0.9, alpha=0.65)
    ax.text(
        0.98,
        0.12,
        rf"$z={z:.1f},\ z'={zbar:.1f},\ Q^2={q2:.0f}\ \mathrm{{GeV}}^2,\ C_{{nn}}={cnn:g}$",
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=10,
    )
    ax.legend(loc="lower left", frameon=True, fontsize=9)
    ax.set_xlim(min(pperp_grid), max(pperp_grid))
    ax.set_title("JAM3D Collins spin correlation", pad=10)
    ax.set_xlabel(r"$P_\perp = P'_\perp\ \mathrm{(GeV)}$")
    ax.set_ylabel(r"$\langle \sin\phi\,\sin\phi'\rangle$")
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
    """Calculate the spin-correlation observable and write CSV/PNG outputs."""
    Z = 0.3
    ZBAR = 0.3
    Q2 = 4.0
    CNN = -1.0
    MAX_KT = 1.0
    PPERP_GRID = linspace(0.0, Z * MAX_KT, 200)
    CHANNELS = ["fav", "unf"]
    REPLICA_SELECTION = 50
    TAG = "JAM3D_2022"
    CSV_OUTPUT = output_path("spin_correlation.csv")
    PLOT_OUTPUT = output_path("spin_correlation.png")

    tmd, replicas, jam3dlib_path, jam3d_dev_lib_path = load_collins_fit(
        tag=TAG,
        replica_selection=REPLICA_SELECTION,
    )
    rows = spin_correlation_rows(
        tmd=tmd,
        pperp_grid=PPERP_GRID,
        z=Z,
        zbar=ZBAR,
        q2=Q2,
        channels=CHANNELS,
        replicas=replicas,
        cnn=CNN,
    )
    write_rows(
        rows=rows,
        output=CSV_OUTPUT,
        fieldnames=[
            "P_perp",
            "Pbar_perp",
            "z",
            "zbar",
            "Q2",
            "Cnn",
            "mh",
            "channel",
            "flavor_index",
            "nrep",
            "coefficient_mean",
            "coefficient_std",
            "spin_correlation_mean",
            "spin_correlation_std",
            "D1_mean",
            "D1_std",
            "H1perp_mean",
            "H1perp_std",
        ],
    )
    plot_rows(rows, PPERP_GRID, CHANNELS, Z, ZBAR, Q2, CNN, PLOT_OUTPUT)

    print(
        f"Wrote {CSV_OUTPUT} and {PLOT_OUTPUT} for z={Z}, zbar={ZBAR}, "
        f"Q2={Q2}, Cnn={CNN}, {len(CHANNELS)} channels, and "
        f"{len(replicas)} replicas.\n"
        f"jam3dlib: {jam3dlib_path}\n"
        f"jam3d_dev_lib: {jam3d_dev_lib_path}"
    )


if __name__ == "__main__":
    main()
