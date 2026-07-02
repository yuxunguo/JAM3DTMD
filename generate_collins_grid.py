"""Generate a Collins TMD CSV grid from JAM3D."""

from __future__ import annotations

from time import perf_counter

from Collins import collins_tmd_grid_rows, linspace, load_collins_fit, output_path
from Collins import write_rows


def format_elapsed(seconds: float) -> str:
    """Format elapsed seconds as a compact human-readable duration."""
    minutes, seconds = divmod(seconds, 60.0)
    if minutes >= 1:
        return f"{int(minutes)}m {seconds:.1f}s"
    return f"{seconds:.1f}s"


def main() -> None:
    """Generate a CSV grid of favored/unfavored pion Collins TMD values."""
    total_start = perf_counter()
    Q2_VALUES = [4.0, 16.0]
    Z_GRID = linspace(0.01, 0.99, 20)
    PT_GRID = linspace(0.01, 0.99, 20)
    FLAVORS = ["fav", "unf"]
    REPLICA_SELECTION = 50
    TAG = "JAM3D_2022"
    OUTPUT = output_path("collins_grid.csv")
    WEIGHTED = False

    load_start = perf_counter()
    tmd, replicas, jam3dlib_path, jam3d_dev_lib_path = load_collins_fit(
        tag=TAG,
        replica_selection=REPLICA_SELECTION,
    )
    print(f"Loaded {TAG} in {format_elapsed(perf_counter() - load_start)}.")

    rows = []
    for q2 in Q2_VALUES:
        q2_start = perf_counter()
        q2_rows = collins_tmd_grid_rows(
            tmd=tmd,
            z_grid=Z_GRID,
            pt_grid=PT_GRID,
            q2=q2,
            flavors=FLAVORS,
            replicas=replicas,
            weighted=WEIGHTED,
        )
        rows.extend(q2_rows)
        print(
            f"Generated Q2={q2:g} grid with {len(q2_rows)} rows in "
            f"{format_elapsed(perf_counter() - q2_start)}."
        )

    write_start = perf_counter()
    output_rows = [
        {
            "x": row["z"],
            "kT": row["pT"],
            "Q2": row["Q2"],
            "flavor": row["flavor"],
            "flavor_index": row["flavor_index"],
            "nrep": row["nrep"],
            "collins_mean": row["mean"],
            "collins_std": row["std"],
        }
        for row in rows
    ]
    write_rows(
        rows=output_rows,
        output=OUTPUT,
        fieldnames=[
            "x",
            "kT",
            "Q2",
            "flavor",
            "flavor_index",
            "nrep",
            "collins_mean",
            "collins_std",
        ],
    )
    write_elapsed = perf_counter() - write_start
    total_elapsed = perf_counter() - total_start

    print(
        f"Wrote {OUTPUT} with {len(Z_GRID)} z-points, "
        f"{len(PT_GRID)} pT-points, {len(FLAVORS)} flavors, "
        f"Q2={Q2_VALUES}, "
        f"and {len(replicas)} replicas.\n"
        f"write time: {format_elapsed(write_elapsed)}\n"
        f"total time: {format_elapsed(total_elapsed)}\n"
        f"jam3dlib: {jam3dlib_path}\n"
        f"jam3d_dev_lib: {jam3d_dev_lib_path}"
    )


if __name__ == "__main__":
    main()
