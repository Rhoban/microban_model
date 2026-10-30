"""Put a cut climbing foot into this model: the drawn mesh, its convex collision pieces, and the
foot's mass properties, on both climbing feet, in robot_climb.xml and robot_climb.urdf.

    uv run --with numpy --with coacd tools/cut_foot.py FOOT.stl [--rotate-x180]

Run it from this repository's root, on the cut-foot branch, with git-lfs installed. FOOT.stl: the
cut foot, a closed binary STL in metres, in the frame of Rhoban's foot__configuration_right.stl
(--rotate-x180 if it is exported turned 180 deg about x, as from Onshape). Then commit, push, and
pin the new commit where the model is used.

Mass properties: Rhoban's CAD values (RHOBAN, below) minus the removed part, i.e. the difference
between Rhoban's mesh and FOOT.stl at the foot's mean density (CAD mass / Rhoban mesh volume), so
whatever the CAD has beyond a uniform solid (screws, infill) stays. Collision pieces: CoACD,
threshold 0.02, at most 64 vertices each, seed 0 (another machine can give slightly other pieces).
"""

import argparse
import re
import struct
import subprocess
from pathlib import Path

import numpy as np

RHOBAN = "8e42ef70673b5747dbdca3efaf60bc0f8ec7ea6d"  # the commit the fork started from
NAME = "foot__configuration_right"
ASSETS = Path("src/microban_model/mjcf/assets")
MJCF = Path("src/microban_model/mjcf/robot_climb.xml")
URDF = Path("src/microban_model/urdf/robot_climb.urdf")
FEET = ("foot_climb", "foot_climb_2")
# Rhoban's inertial of each climbing foot, body frame (robot_climb.xml at RHOBAN)
MASS = 0.0299799
COM = np.array([-1.0695e-08, -0.0115199, 0.0605776])
QUAT = (0.705504, -0.0475801, 0.0475801, 0.705504)
DIAG = np.array([4.07283e-05, 4.01815e-05, 5.73401e-06])
# where the foot body places the foot mesh (the geoms' pos and quat in robot_climb.xml)
GEOM_POS = np.array([0.0, -0.0095, 0.0385])
GEOM_QUAT = (0.0, 0.0, 0.707107, 0.707107)


def read_stl(data: bytes) -> tuple[np.ndarray, np.ndarray]:
    (n,) = struct.unpack("<I", data[80:84])
    if len(data) != 84 + 50 * n:
        raise ValueError("not a binary STL")
    rec = np.dtype([("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    corners = np.frombuffer(data, dtype=rec, count=n, offset=84)["v"].reshape(-1, 3)
    keys, index = np.unique(np.round(corners.astype(float) * 1e7).astype(np.int64), axis=0, return_inverse=True)
    return keys * 1e-7, index.reshape(-1, 3)


def write_stl(path: Path, v: np.ndarray, t: np.ndarray, header: bytes) -> None:
    tri = v[t]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-30
    rec = np.zeros(len(t), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    rec["n"], rec["v"] = n, tri
    path.write_bytes(header.ljust(80, b" ")[:80] + struct.pack("<I", len(t)) + rec.tobytes())


def closed(t: np.ndarray) -> bool:
    edges = np.sort(np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]), axis=1)
    return bool((np.unique(edges, axis=0, return_counts=True)[1] == 2).all())


def quat_matrix(w, x, y, z):
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])  # fmt: skip


def body_props(v, t):
    """(volume, centroid, integral of x x^T dV), uniform solid, in the foot body's frame."""
    v = np.einsum("ij,nj->ni", quat_matrix(*GEOM_QUAT), v) + GEOM_POS
    a, b, c = v[t[:, 0]], v[t[:, 1]], v[t[:, 2]]
    vol6 = np.einsum("ij,ij->i", a, np.cross(b, c))
    volume = vol6.sum() / 6
    centroid = ((a + b + c) * vol6[:, None]).sum(0) / 24 / volume
    s = a + b + c
    second = np.einsum("n,ni,nj->ij", vol6 / 120, s, s)
    second += sum(np.einsum("n,ni,nj->ij", vol6 / 120, p, p) for p in (a, b, c))
    return volume, centroid, second


def inertial(cut_v, cut_t, rhoban_v, rhoban_t):
    """(mass, centre of mass, principal axes, principal moments) of the cut foot, body frame."""
    v0, c0, s0 = body_props(rhoban_v, rhoban_t)
    v1, c1, s1 = body_props(cut_v, cut_t)
    rho = MASS / v0
    removed = rho * (v0 - v1)
    i_cad = quat_matrix(*QUAT) @ np.diag(DIAG) @ quat_matrix(*QUAT).T
    second = (np.trace(i_cad) / 2 * np.eye(3) - i_cad) + MASS * np.outer(COM, COM)
    mass = MASS - removed
    com = (MASS * COM - rho * (v0 * c0 - v1 * c1)) / mass
    s = second - rho * (s0 - s1) - mass * np.outer(com, com)
    moments, axes = np.linalg.eigh(np.trace(s) * np.eye(3) - s)
    order = np.argsort(-moments)
    moments, axes = moments[order], axes[:, order]
    if np.linalg.det(axes) < 0:
        axes[:, 2] *= -1
    cad = quat_matrix(*QUAT)  # of the right-handed sign choices, the one closest to the CAD's
    flips = [np.diag(f) for f in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))]
    axes = max((axes @ f for f in flips), key=lambda r: np.trace(r.T @ cad))
    return mass, com, axes, moments


def quat_of(r):
    w = np.sqrt(max(0.0, 1 + r[0, 0] + r[1, 1] + r[2, 2])) / 2
    x = np.copysign(np.sqrt(max(0.0, 1 + r[0, 0] - r[1, 1] - r[2, 2])) / 2, r[2, 1] - r[1, 2])
    y = np.copysign(np.sqrt(max(0.0, 1 - r[0, 0] + r[1, 1] - r[2, 2])) / 2, r[0, 2] - r[2, 0])
    z = np.copysign(np.sqrt(max(0.0, 1 - r[0, 0] - r[1, 1] + r[2, 2])) / 2, r[1, 0] - r[0, 1])
    q = np.array([w, x, y, z])
    return q / np.linalg.norm(q)


def rpy_of(r):  # URDF: r = Rz(yaw) Ry(pitch) Rx(roll)
    return np.array([np.arctan2(r[2, 1], r[2, 2]), np.arcsin(-np.clip(r[2, 0], -1, 1)), np.arctan2(r[1, 0], r[0, 0])])


def g(x: float) -> str:
    return f"{x:.6g}"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("foot", type=Path)
    parser.add_argument("--rotate-x180", action="store_true", help="the STL is turned 180 deg about x")
    args = parser.parse_args()
    import coacd

    blob = subprocess.run(["git", "show", f"{RHOBAN}:{ASSETS / NAME}.stl"], check=True, capture_output=True).stdout
    rhoban = subprocess.run(["git", "lfs", "smudge"], input=blob, check=True, capture_output=True).stdout
    rhoban_v, rhoban_t = read_stl(rhoban)
    v, t = read_stl(args.foot.read_bytes())
    if args.rotate_x180:
        v = v * [1, -1, -1]
    if not closed(t):
        raise SystemExit(f"{args.foot} is not one closed surface: every edge must join two triangles")
    gap = np.abs(np.r_[v.min(0), v.max(0)] - np.r_[rhoban_v.min(0), rhoban_v.max(0)]).max()
    if gap > 0.005:  # the cut changes some bounds, but not all of them by much
        print(f"warning: its bounds differ from Rhoban's foot by up to {gap * 1e3:.1f} mm: check the frame")

    # written first, then everything from the file as written (float32): the same pieces each run
    write_stl(ASSETS / f"{NAME}.stl", v, t, b"climbing foot, bottom plate's end cut off (rock-climbing-workshop)")
    v, t = read_stl((ASSETS / f"{NAME}.stl").read_bytes())
    mass, com, axes, moments = inertial(v, t, rhoban_v, rhoban_t)
    q, rpy = quat_of(axes), rpy_of(axes)

    for old in ASSETS.glob(f"{NAME}_[0-9][0-9].stl"):
        old.unlink()
    coacd.set_log_level("error")
    parts = coacd.run_coacd(coacd.Mesh(v, t), threshold=0.02, seed=0, decimate=True, max_ch_vertex=64)
    for i, (pv, pt) in enumerate(parts):
        write_stl(ASSETS / f"{NAME}_{i:02d}.stl", np.asarray(pv, float), np.asarray(pt), b"convex piece (CoACD)")

    text = MJCF.read_text()
    meshes = re.search(rf'(^( *)<mesh name="coll_{NAME}_\d\d"[^\n]*\n)+', text, re.M)
    indent = meshes.group(2)
    block = "".join(
        f'{indent}<mesh name="coll_{NAME}_{i:02d}" content_type="model/stl" file="{NAME}_{i:02d}.stl" />\n'
        for i in range(len(parts))
    )
    text = text[: meshes.start()] + block + text[meshes.end() :]

    def geoms(match):
        line = match.group(0).splitlines(keepends=True)[0]
        return "".join(re.sub(rf'mesh="coll_{NAME}_\d\d"', f'mesh="coll_{NAME}_{i:02d}"', line) for i in range(len(parts)))

    text, n = re.subn(rf'(^ *<geom [^\n]*mesh="coll_{NAME}_\d\d" />\n)+', geoms, text, flags=re.M)
    assert n == len(FEET), n
    new = (
        f'<inertial pos="{g(com[0])} {g(com[1])} {g(com[2])}" quat="{g(q[0])} {g(q[1])} {g(q[2])} {g(q[3])}" '
        f'mass="{g(mass)}" diaginertia="{g(moments[0])} {g(moments[1])} {g(moments[2])}" />'
    )
    for foot in FEET:
        text, n = re.subn(rf'(<body name="{foot}"[^\n]*\n\s*)<inertial [^\n]*/>', rf"\g<1>{new}", text)
        assert n == 1, foot
    MJCF.write_text(text)

    urdf = URDF.read_text()
    new = (
        f'<origin xyz="{g(com[0])} {g(com[1])} {g(com[2])}" rpy="{rpy[0]:.8f} {rpy[1]:.8f} {rpy[2]:.8f}"/>\n'
        f'      <mass value="{g(mass)}"/>\n'
        f'      <inertia ixx="{g(moments[0])}" ixy="0" ixz="0" iyy="{g(moments[1])}" iyz="0" izz="{g(moments[2])}"/>'
    )
    for foot in FEET:
        pattern = rf'(<link name="{foot}">\s*<inertial>\s*)<origin [^\n]*/>\s*<mass [^\n]*/>\s*<inertia [^\n]*/>'
        urdf, n = re.subn(pattern, lambda m: m.group(1) + new, urdf)
        assert n == 1, foot
    URDF.write_text(urdf)
    print(f"foot: {MASS * 1e3:.2f} g -> {mass * 1e3:.2f} g; {len(parts)} collision pieces; written")


if __name__ == "__main__":
    main()
