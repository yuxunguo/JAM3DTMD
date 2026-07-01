"""Generate a Collins-function grid from the JAM3D TMD library.

The notebook export in ``copy_of_jam3d_library.py`` evaluates the pion
Collins TMD as:

    tmd.eval(z, Q2, kT, "pi", "collinspi", irep, icol=False)

This script turns that call into a reusable CSV grid in x/z and kT for a fixed
Q2 and selected flavors.
"""

from __future__ import annotations

import csv
import importlib
import math
import os
import sys
import types
from pathlib import Path


FLAVOR_INDEX = {
    "g": 0,
    "u": 1,
    "ub": 2,
    "d": 3,
    "db": 4,
    "s": 5,
    "sb": 6,
    "c": 7,
    "cb": 8,
    "b": 9,
    "bb": 10,
    "fav": 1,
    "unf": 3,
}

MPI_GEV = 0.135
REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_JAM3DLIB_PATHS = (
    REPO_ROOT / "jam3dlib",
    REPO_ROOT / "jam3dlib-master",
)
DEFAULT_JAM3D_DEV_LIB_PATHS = (
    REPO_ROOT / "jam3d_dev_lib",
    REPO_ROOT / "jam3d_dev_lib-main",
)


def linspace(start: float, stop: float, count: int) -> list[float]:
    if count < 1:
        raise ValueError("grid count must be positive")
    if count == 1:
        return [start]
    step = (stop - start) / (count - 1)
    return [start + i * step for i in range(count)]


def select_replicas(
    selection: int | tuple[int, int] | list[int] | str,
    nrep: int,
) -> list[int]:
    """Select replicas as all, N, (start, stop), or explicit index list."""
    if selection == "all":
        return list(range(nrep))
    if isinstance(selection, int):
        return list(range(min(selection, nrep)))
    if isinstance(selection, tuple):
        start, stop = selection
        return list(range(start, min(stop, nrep)))
    return [irep for irep in selection if 0 <= irep < nrep]


def validate_flavors(flavors: list[str]) -> list[str]:
    unknown = [flavor for flavor in flavors if flavor not in FLAVOR_INDEX]
    if unknown:
        names = ", ".join(sorted(FLAVOR_INDEX))
        raise ValueError(
            f"unknown flavor(s): {', '.join(unknown)}. Choices: {names}"
        )
    return flavors


def path_from_env(name: str) -> Path | None:
    value = os.environ.get(name)
    return Path(value).expanduser() if value else None


def resolve_library_path(
    label: str,
    configured: Path | None,
    env_name: str,
    candidates: tuple[Path, ...],
    marker_file: str,
) -> Path:
    paths = [path for path in (path_from_env(env_name), configured) if path is not None]
    paths.extend(candidates)

    checked: list[Path] = []
    for path in paths:
        path = path.resolve()
        checked.append(path)
        if (path / marker_file).exists():
            return path

    checked_paths = "\n".join(f"  {path}" for path in checked)
    raise SystemExit(
        f"Could not find {label}.\n"
        f"Checked:\n{checked_paths}\n\n"
        f"Set {env_name} or update the path constants in main()."
    )


def install_jam3d_compatibility_shims() -> None:
    from scipy import special

    class _ZmqUnavailableError(RuntimeError):
        pass

    class _ZmqContextFallback:
        def socket(self, *_args, **_kwargs):
            raise _ZmqUnavailableError(
                "pyzmq is not installed; JAM3D parallel workers are unavailable."
            )

    class _MpmathFpFallback:
        euler = 0.5772156649015329

        @staticmethod
        def psi(order, value):
            if order == 0:
                return special.digamma(value)
            return special.polygamma(order, value)

        @staticmethod
        def zeta(value):
            if value == 2:
                return math.pi**2 / 6.0
            if value == 3:
                return 1.2020569031595942
            return special.zeta(value, 1.0)

    zmq = types.ModuleType("zmq")
    zmq.REP = object()
    zmq.REQ = object()
    zmq.ZMQError = _ZmqUnavailableError
    zmq.Context = _ZmqContextFallback
    sys.modules["zmq"] = zmq

    mpmath = types.ModuleType("mpmath")
    mpmath.fp = _MpmathFpFallback()
    mpmath.mp = mpmath.fp
    sys.modules["mpmath"] = mpmath

    aux = importlib.import_module("qcdlib._aux")
    qcdlib = importlib.import_module("qcdlib")
    sys.modules["qcdlib.aux"] = aux
    setattr(qcdlib, "aux", aux)


def import_tmd(jam3dlib: Path | None, jam3d_dev_lib: Path | None):
    jam3dlib = resolve_library_path(
        label="jam3dlib",
        configured=jam3dlib,
        env_name="JAM3DLIB_PATH",
        candidates=DEFAULT_JAM3DLIB_PATHS,
        marker_file="tmd.py",
    )
    jam3d_dev_lib = resolve_library_path(
        label="jam3d_dev_lib",
        configured=jam3d_dev_lib,
        env_name="JAM3D_DEV_LIB_PATH",
        candidates=DEFAULT_JAM3D_DEV_LIB_PATHS,
        marker_file="tools/config.py",
    )

    os.environ["JAM3D"] = str(jam3d_dev_lib)
    for path in (jam3d_dev_lib, jam3dlib):
        path_text = str(path)
        if path_text not in sys.path:
            sys.path.insert(0, path_text)
    install_jam3d_compatibility_shims()

    from tmd import TMD
    return TMD, jam3dlib, jam3d_dev_lib


def load_tmd(TMD, tag: str, jam3dlib_path: Path):
    old_cwd = Path.cwd()
    os.chdir(jam3dlib_path)
    tmd = TMD(tag)
    os.chdir(old_cwd)
    return tmd


def build_rows(
    tmd,
    x_grid: list[float],
    kt_grid: list[float],
    q2: float,
    flavors: list[str],
    replicas: list[int],
    weighted: bool,
):
    from tools.config import conf

    collins = conf["collinspi"]
    points = [(x, kt) for x in x_grid for kt in kt_grid]
    flavor_indices = [(flavor, FLAVOR_INDEX[flavor]) for flavor in flavors]
    sums = [[0.0 for _ in flavor_indices] for _ in points]
    sums_sq = [[0.0 for _ in flavor_indices] for _ in points]

    for irep in replicas:
        tmd.parman.set_new_params(tmd.par[irep], initial=True)
        for point_index, (x, kt) in enumerate(points):
            values = collins.get_tmd(x, q2, kt, "pi", "collinspi", icol=False)
            weight = x**2 * kt**2 / (2.0 * MPI_GEV**2) if weighted else 1.0
            for flavor_position, (_flavor, flavor_index) in enumerate(flavor_indices):
                value = float(values[flavor_index]) * weight
                sums[point_index][flavor_position] += value
                sums_sq[point_index][flavor_position] += value * value

    nrep = len(replicas)
    for point_index, (x, kt) in enumerate(points):
        for flavor_position, (flavor, flavor_index) in enumerate(flavor_indices):
            mean = sums[point_index][flavor_position] / nrep
            mean_sq = sums_sq[point_index][flavor_position] / nrep
            variance = max(mean_sq - mean * mean, 0.0)
            yield {
                "x": x,
                "kT": kt,
                "Q2": q2,
                "flavor": flavor,
                "flavor_index": flavor_index,
                "nrep": nrep,
                "collins_mean": mean,
                "collins_std": math.sqrt(variance),
            }


def main() -> None:
    # Change these variables for a new grid.
    Q2 = 4.0
    X_GRID = linspace(0.01, 0.99, 100)
    KT_GRID = linspace(0.01, 0.99, 100)
    FLAVORS = ["fav", "unf"]
    REPLICA_SELECTION = 50  # use "all", an integer N, (start, stop), or [indices]
    TAG = "JAM3D_2022"
    OUTPUT = Path("collins_grid.csv")
    JAM3DLIB_PATH = None
    JAM3D_DEV_LIB_PATH = None
    WEIGHTED = False  # True saves x^2 kT^2/(2 Mpi^2) times the raw Collins TMD

    flavors = validate_flavors(FLAVORS)
    TMD, jam3dlib_path, jam3d_dev_lib_path = import_tmd(
        JAM3DLIB_PATH,
        JAM3D_DEV_LIB_PATH,
    )
    tmd = load_tmd(TMD, TAG, jam3dlib_path)
    replicas = select_replicas(REPLICA_SELECTION, tmd.nrep)
    if not replicas:
        raise SystemExit("No valid replicas selected.")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "x",
        "kT",
        "Q2",
        "flavor",
        "flavor_index",
        "nrep",
        "collins_mean",
        "collins_std",
    ]
    with OUTPUT.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(
            build_rows(
                tmd=tmd,
                x_grid=X_GRID,
                kt_grid=KT_GRID,
                q2=Q2,
                flavors=flavors,
                replicas=replicas,
                weighted=WEIGHTED,
            )
        )

    print(
        f"Wrote {OUTPUT} with {len(X_GRID)} x-points, "
        f"{len(KT_GRID)} kT-points, {len(flavors)} flavors, "
        f"and {len(replicas)} replicas.\n"
        f"jam3dlib: {jam3dlib_path}\n"
        f"jam3d_dev_lib: {jam3d_dev_lib_path}"
    )


if __name__ == "__main__":
    main()
