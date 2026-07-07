"""Reusable JAM3D Collins-function utilities."""

from __future__ import annotations

import csv
import importlib
import math
import os
import sys
import types
from pathlib import Path
from typing import Callable, Sequence


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
OUTPUT_DIR = REPO_ROOT / "Output"
DEFAULT_JAM3DLIB_PATHS = (
    REPO_ROOT / "jam3dlib",
    REPO_ROOT / "jam3dlib-master",
)
DEFAULT_JAM3D_DEV_LIB_PATHS = (
    REPO_ROOT / "jam3d_dev_lib",
    REPO_ROOT / "jam3d_dev_lib-main",
)


def linspace(start: float, stop: float, count: int) -> list[float]:
    """Return count evenly spaced points including both endpoints."""
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
    """Expand a replica selection into valid JAM3D replica indices."""
    if selection == "all":
        return list(range(nrep))
    if isinstance(selection, int):
        return list(range(min(selection, nrep)))
    if isinstance(selection, tuple):
        start, stop = selection
        return list(range(start, min(stop, nrep)))
    return [irep for irep in selection if 0 <= irep < nrep]


def validate_flavors(flavors: list[str]) -> list[str]:
    """Validate flavor labels used by this analysis and return them unchanged."""
    unknown = [flavor for flavor in flavors if flavor not in FLAVOR_INDEX]
    if unknown:
        names = ", ".join(sorted(FLAVOR_INDEX))
        raise ValueError(
            f"unknown flavor(s): {', '.join(unknown)}. Choices: {names}"
        )
    return flavors


def write_rows(
    rows: list[dict[str, float | int | str]],
    output: Path,
    fieldnames: list[str],
) -> None:
    """Write dictionaries to a CSV file, creating the output directory first."""
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def output_path(filename: str) -> Path:
    """Return a path inside the repository Output directory."""
    return OUTPUT_DIR / filename


def path_from_env(name: str) -> Path | None:
    """Read a filesystem path from an environment variable if it is set."""
    value = os.environ.get(name)
    return Path(value).expanduser() if value else None


def resolve_library_path(
    label: str,
    configured: Path | None,
    env_name: str,
    candidates: tuple[Path, ...],
    marker_file: str,
) -> Path:
    """Resolve a JAM3D library path from env, explicit config, or local defaults."""
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
    """Install small compatibility modules needed by the bundled JAM3D code."""
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
    """Load JAM3D paths, install shims, and import the public TMD wrapper."""
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
    """Instantiate a JAM3D TMD object from the jam3dlib data directory."""
    old_cwd = Path.cwd()
    os.chdir(jam3dlib_path)
    tmd = TMD(tag)
    os.chdir(old_cwd)
    return tmd


def load_collins_fit(
    tag: str,
    replica_selection: int | tuple[int, int] | list[int] | str,
    jam3dlib_path: Path | None = None,
    jam3d_dev_lib_path: Path | None = None,
):
    """Load a JAM3D fit and return its TMD wrapper plus selected replicas."""
    TMD, resolved_jam3dlib, resolved_jam3d_dev_lib = import_tmd(
        jam3dlib_path,
        jam3d_dev_lib_path,
    )
    tmd = load_tmd(TMD, tag, resolved_jam3dlib)
    replicas = select_replicas(replica_selection, tmd.nrep)
    if not replicas:
        raise SystemExit("No valid replicas selected.")
    return tmd, replicas, resolved_jam3dlib, resolved_jam3d_dev_lib


def flavor_indices(flavors: list[str]) -> list[tuple[str, int]]:
    """Map flavor labels to JAM3D flavor-array indices."""
    return [(flavor, FLAVOR_INDEX[flavor]) for flavor in validate_flavors(flavors)]


def mean_std_rows(
    coordinates: list[dict[str, float]],
    flavors: list[str],
    replicas: list[int],
    value_func: Callable[[dict[str, float], int], Sequence[float]],
) -> list[dict[str, float | int | str]]:
    """Evaluate replica-dependent flavor vectors and summarize mean/std rows."""
    indexed_flavors = flavor_indices(flavors)
    sums = [[0.0 for _ in indexed_flavors] for _ in coordinates]
    sums_sq = [[0.0 for _ in indexed_flavors] for _ in coordinates]

    for irep in replicas:
        for coord_index, coord in enumerate(coordinates):
            values = value_func(coord, irep)
            for flavor_position, (_flavor, flavor_index) in enumerate(indexed_flavors):
                value = float(values[flavor_index])
                sums[coord_index][flavor_position] += value
                sums_sq[coord_index][flavor_position] += value * value

    nrep = len(replicas)
    rows: list[dict[str, float | int | str]] = []
    for coord_index, coord in enumerate(coordinates):
        for flavor_position, (flavor, flavor_index) in enumerate(indexed_flavors):
            mean = sums[coord_index][flavor_position] / nrep
            mean_sq = sums_sq[coord_index][flavor_position] / nrep
            variance = max(mean_sq - mean * mean, 0.0)
            row: dict[str, float | int | str] = dict(coord)
            row.update(
                {
                    "flavor": flavor,
                    "flavor_index": flavor_index,
                    "nrep": nrep,
                    "mean": mean,
                    "std": math.sqrt(variance),
                }
            )
            rows.append(row)
    return rows


def collins_tmd_grid_rows(
    tmd,
    z_grid: list[float],
    pt_grid: list[float],
    q2: float,
    flavors: list[str],
    replicas: list[int],
    weighted: bool,
) -> list[dict[str, float | int | str]]:
    """Build mean/std rows for the pion Collins TMD over a z and pT grid.

    The TMD values are evaluated through JAM3D's public tmd.eval wrapper, using
    the same argument order shown in JAM3D_Library.ipynb.
    """
    coordinates = [{"z": z, "pT": pt, "Q2": q2} for z in z_grid for pt in pt_grid]

    def value_func(coord: dict[str, float], irep: int) -> Sequence[float]:
        z = coord["z"]
        pt = coord["pT"]
        # Explicit JAM3D wrapper call: icol=False returns the full pT-dependent TMD.
        values = tmd.eval(z, q2, pt, "pi", "collinspi", irep, icol=False)
        weight = z**2 * pt**2 / (2.0 * MPI_GEV**2) if weighted else 1.0
        return weight * values

    return mean_std_rows(coordinates, flavors, replicas, value_func)


def z_collinear_collins_rows(
    tmd,
    z_grid: list[float],
    q2: float,
    flavors: list[str],
    replicas: list[int],
) -> list[dict[str, float | int | str]]:
    """Build mean/std rows for z times the collinear pion Collins function."""
    coordinates = [{"z": z, "Q2": q2} for z in z_grid]

    def value_func(coord: dict[str, float], irep: int) -> Sequence[float]:
        z = coord["z"]
        # Explicit note: icol=True returns the collinear get_C(z, Q2) piece.
        # The kT argument is required by tmd.eval but ignored in this mode.
        values = tmd.eval(z, q2, 0.0, "pi", "collinspi", irep, icol=True)
        return z * values

    return mean_std_rows(coordinates, flavors, replicas, value_func)


def fixed_z_weighted_collins_rows(
    tmd,
    z: float,
    pt_grid: list[float],
    q2: float,
    flavors: list[str],
    replicas: list[int],
) -> list[dict[str, float | int | str]]:
    """Build mean/std rows for weighted pion Collins TMD values at fixed z."""
    coordinates = [{"z": z, "pT": pt, "z_pT": z * pt, "Q2": q2} for pt in pt_grid]

    def value_func(coord: dict[str, float], irep: int) -> Sequence[float]:
        pt = coord["pT"]
        # Explicit JAM3D wrapper call: icol=False returns the full pT-dependent TMD.
        values = tmd.eval(z, q2, pt, "pi", "collinspi", irep, icol=False)
        weight = z**2 * pt**2 / (2.0 * MPI_GEV**2)
        return weight * values

    return mean_std_rows(coordinates, flavors, replicas, value_func)


def fixed_z_unpolarized_ff_rows(
    tmd,
    z: float,
    pt_grid: list[float],
    q2: float,
    flavors: list[str],
    replicas: list[int],
) -> list[dict[str, float | int | str]]:
    """Build mean/std rows for unpolarized pion TMD FF values at fixed z."""
    coordinates = [{"z": z, "pT": pt, "z_pT": z * pt, "Q2": q2} for pt in pt_grid]

    def value_func(coord: dict[str, float], irep: int) -> Sequence[float]:
        pt = coord["pT"]
        values = tmd.eval(z, q2, pt, "pi", "ff", irep, icol=False)
        return values

    return mean_std_rows(coordinates, flavors, replicas, value_func)
