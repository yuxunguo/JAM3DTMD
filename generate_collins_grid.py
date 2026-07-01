"""Generate a Collins TMD CSV grid from JAM3D."""

from __future__ import annotations

from Collins import collins_tmd_grid_rows, linspace, load_collins_fit, output_path
from Collins import write_rows


def main() -> None:
    Q2 = 4.0
    Z_GRID = linspace(0.01, 0.99, 100)
    PT_GRID = linspace(0.01, 0.99, 100)
    FLAVORS = ["fav", "unf"]
    REPLICA_SELECTION = 50
    TAG = "JAM3D_2022"
    OUTPUT = output_path("collins_grid.csv")
    WEIGHTED = False

    tmd, replicas, jam3dlib_path, jam3d_dev_lib_path = load_collins_fit(
        tag=TAG,
        replica_selection=REPLICA_SELECTION,
    )
    rows = collins_tmd_grid_rows(
        tmd=tmd,
        z_grid=Z_GRID,
        pt_grid=PT_GRID,
        q2=Q2,
        flavors=FLAVORS,
        replicas=replicas,
        weighted=WEIGHTED,
    )
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

    print(
        f"Wrote {OUTPUT} with {len(Z_GRID)} z-points, "
        f"{len(PT_GRID)} pT-points, {len(FLAVORS)} flavors, "
        f"and {len(replicas)} replicas.\n"
        f"jam3dlib: {jam3dlib_path}\n"
        f"jam3d_dev_lib: {jam3d_dev_lib_path}"
    )


if __name__ == "__main__":
    main()
