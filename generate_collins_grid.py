"""Generate a Collins TMD CSV grid from JAM3D."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from math import ceil
from os import cpu_count
from time import perf_counter

from Collins import collins_tmd_grid_rows, linspace, load_collins_fit, output_path
from Collins import write_rows


_WORKER_TMD = None
_WORKER_REPLICAS: list[int] = []
_WORKER_JAM3DLIB_PATH = None
_WORKER_JAM3D_DEV_LIB_PATH = None


def format_elapsed(seconds: float) -> str:
    """Format elapsed seconds as a compact human-readable duration."""
    minutes, seconds = divmod(seconds, 60.0)
    if minutes >= 1:
        return f"{int(minutes)}m {seconds:.1f}s"
    return f"{seconds:.1f}s"


def chunked(values: list[float], nchunks: int) -> list[list[float]]:
    """Split values into at most nchunks non-empty chunks."""
    chunk_size = max(1, ceil(len(values) / nchunks))
    return [values[i : i + chunk_size] for i in range(0, len(values), chunk_size)]


def init_worker(tag: str, replica_selection: int | tuple[int, int] | list[int] | str):
    """Load the JAM3D fit once inside each worker process."""
    global _WORKER_TMD
    global _WORKER_REPLICAS
    global _WORKER_JAM3DLIB_PATH
    global _WORKER_JAM3D_DEV_LIB_PATH

    (
        _WORKER_TMD,
        _WORKER_REPLICAS,
        _WORKER_JAM3DLIB_PATH,
        _WORKER_JAM3D_DEV_LIB_PATH,
    ) = load_collins_fit(tag=tag, replica_selection=replica_selection)


def generate_q2_z_chunk(
    task: tuple[float, list[float], list[float], list[str], bool],
) -> dict[str, object]:
    """Generate one Q2/z-grid chunk using the worker-local JAM3D TMD object."""
    if _WORKER_TMD is None:
        raise RuntimeError("worker TMD fit was not initialized")

    q2, z_grid, pt_grid, flavors, weighted = task
    start = perf_counter()
    rows = collins_tmd_grid_rows(
        tmd=_WORKER_TMD,
        z_grid=z_grid,
        pt_grid=pt_grid,
        q2=q2,
        flavors=flavors,
        replicas=_WORKER_REPLICAS,
        weighted=weighted,
    )
    return {
        "q2": q2,
        "rows": rows,
        "z_count": len(z_grid),
        "row_count": len(rows),
        "elapsed": perf_counter() - start,
        "nrep": len(_WORKER_REPLICAS),
        "jam3dlib_path": _WORKER_JAM3DLIB_PATH,
        "jam3d_dev_lib_path": _WORKER_JAM3D_DEV_LIB_PATH,
    }


def main() -> None:
    """Generate a CSV grid of favored/unfavored pion Collins TMD values."""
    total_start = perf_counter()
    Q2_VALUES = [4.0, 16.0]
    Z_GRID = linspace(0.01, 0.99, 50)
    PT_GRID = linspace(0.01, 0.99, 50)
    FLAVORS = ["fav", "unf"]
    REPLICA_SELECTION = 50
    TAG = "JAM3D_2022"
    OUTPUT = output_path("collins_grid.csv")
    WEIGHTED = False
    MAX_WORKERS = min(4, cpu_count() or 1)

    workers = max(1, min(MAX_WORKERS, len(Q2_VALUES) * len(Z_GRID)))
    z_chunks = chunked(Z_GRID, workers)
    tasks = [
        (q2, z_chunk, PT_GRID, FLAVORS, WEIGHTED)
        for q2 in Q2_VALUES
        for z_chunk in z_chunks
    ]
    print(
        f"Generating {len(tasks)} chunks with {workers} worker processes "
        f"for Q2={Q2_VALUES}."
    )

    chunk_results = []
    rows = []
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=init_worker,
        initargs=(TAG, REPLICA_SELECTION),
    ) as executor:
        future_to_task = {
            executor.submit(generate_q2_z_chunk, task): task for task in tasks
        }
        for future in as_completed(future_to_task):
            q2, z_chunk, _pt_grid, _flavors, _weighted = future_to_task[future]
            result = future.result()
            chunk_results.append(result)
            rows.extend(result["rows"])
            print(
                f"Generated Q2={q2:g}, z chunk {z_chunk[0]:.3f}-{z_chunk[-1]:.3f} "
                f"with {result['row_count']} rows in "
                f"{format_elapsed(float(result['elapsed']))}."
            )

    rows.sort(key=lambda row: (float(row["Q2"]), float(row["z"]), float(row["pT"])))
    q2_elapsed = {
        q2: sum(
            float(result["elapsed"])
            for result in chunk_results
            if float(result["q2"]) == q2
        )
        for q2 in Q2_VALUES
    }
    nrep = int(chunk_results[0]["nrep"]) if chunk_results else 0
    jam3dlib_path = chunk_results[0]["jam3dlib_path"] if chunk_results else ""
    jam3d_dev_lib_path = (
        chunk_results[0]["jam3d_dev_lib_path"] if chunk_results else ""
    )
    for q2 in Q2_VALUES:
        q2_rows = [row for row in rows if float(row["Q2"]) == q2]
        print(
            f"Finished Q2={q2:g} with {len(q2_rows)} rows; "
            f"worker CPU time {format_elapsed(q2_elapsed[q2])}."
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
        f"and {nrep} replicas.\n"
        f"write time: {format_elapsed(write_elapsed)}\n"
        f"total time: {format_elapsed(total_elapsed)}\n"
        f"jam3dlib: {jam3dlib_path}\n"
        f"jam3d_dev_lib: {jam3d_dev_lib_path}"
    )


if __name__ == "__main__":
    main()
