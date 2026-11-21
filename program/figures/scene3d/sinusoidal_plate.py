#!/usr/bin/env python3
"""Problem figure of the simply supported square plate under a sinusoidal load.

The figure is a 2.5D isometric sketch built with ``mechanicskit.scene3d``: the
plate, the eight hinge supports around its edges (an ear plate under the edge,
pinned through the cheeks of a clevis that stands on a shoe) and the load arrows
are solids and true 3D annotations, projected by one orthographic camera.  Everything follows from the
parameters of Section "Simply support rectangular plate with sinusoidal
transverse load": a square plate of side ``L``, simply supported on all four
edges, loaded by ``q = sin(pi x) sin(pi y)``.

Beyond the far edges of the plate stand its two 1D projection models: the strip
seen along ``y`` is a simply supported beam under ``sin(pi x)``, drawn in the
plane ``y = L + d1``, and the twin along ``x`` carries ``sin(pi y)`` in the plane
``x = L + d1``.  Each is drawn as a structural calculation diagram, only tilted
with the plate like the rest of the scene: the beam is one line on its
mid-surface, the load a row of arrows under its sine envelope, and the supports
are line work -- a small hinge circle at each beam end, a triangle down to the
hatched ground at the pinned end and the same triangle on two rollers at the
other end.  The hinge points are white discs, so the beam line is not seen
inside them, and the ground is a filled bar rather than a hatched line.

    python3 program/figures/scene3d/sinusoidal_plate.py

Output:
    manuscript/figs/sinusoidal_plate/{sinusoidal_plate.json, report.json,
        az-35_el32.svg/.png, ...}
    manuscript/figs/sinusoidal_plate.png, .pdf   the figure of Section 5.1

The camera is az = -35 deg, el = 32 deg: an axonometric view from the near
corner, x running to the lower right, y to the lower left, so both dimension
lines stay visible.

Dependencies (the venv of the manuscript has both):

    uv pip install "mechanicskit[cad]"          # build123d + OpenCascade
"""
from __future__ import annotations

import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np

from mechanicskit import scene3d as s3
from mechanicskit import sketch as sk

# --------------------------------------------------------------- parameters
L = 1000.0        # plate side, mm (the model is L x L, the paper uses (0,1)^2)
t = 30.0          # plate thickness, mm
dis_t = 8.0       # thickness of the white hinge discs of the 1D models, mm
gnd_h, gnd_t = 18.0, 6.0    # height and thickness of the filled ground bar, mm
n = 6             # load arrows per side
l0, l1 = 70.0, 190.0   # load arrow, mm: l0 = shortest, so every tail shows; l1 carries q

d1 = 260.0        # distance of each 1D model outside the far plate edge, mm
m = 9             # load arrows of a 1D model
b0, b1 = 35.0, 100.0   # its load arrow, mm

# the line work of the 1D support symbols, as in a structural calculation diagram
g1 = -73.0        # ground line, mm
r_h = 14.0        # hinge circle at the beam ends, mm
r_r = 7.0         # roller circles of the movable end, mm
w_tri = 110.0     # base of the support triangles, and the width of the ground bar
z_apex = t/2 - r_h - 4.0    # the triangle apex just clears the hinge disc
z_base_pin = g1             # the pinned triangle stands on the ground line itself
z_base_rol = g1 + 2*r_r     # the roller triangle stands on its two rollers

ROOT = Path(__file__).resolve().parents[3]
FIGDIR = ROOT / "manuscript" / "figs"
OUT = FIGDIR / "sinusoidal_plate"

# ------------------------------------------------------------------- solids
# The plate spans the domain (0, L) x (0, L) at z in (0, t).  It carries no
# support of its own: the boundary conditions are read off the two 1D models,
# whose supports are drawn there.  So the only solids are the plate, and for the
# 1D models a white disc at each hinge (the beam line stops at it, as a hinge
# symbol is drawn) and the filled ground bar its triangle stands on.
parts = [
    {"id": "plate", "kind": "plate", "color": "table",
     "args": {"lo": [0, 0, 0], "hi": [L, L, t]}},
]


def box_along(pid, c, u, n, size_u, size_n, z0, z1, color="steel"):
    """An axis-aligned box at c, size_u along u and size_n along n, from z0 to z1."""
    c, u, n = (np.asarray(v, float) for v in (c, u, n))
    h = np.vstack([c + s1*size_u/2*u + s2*size_n/2*n
                   for s1 in (-1, 1) for s2 in (-1, 1)])
    parts.append({"id": pid, "kind": "box", "color": color,
                  "args": {"lo": [h[:, 0].min(), h[:, 1].min(), z0],
                           "hi": [h[:, 0].max(), h[:, 1].max(), z1]}})


# -------------------------------------------------------------- annotations
ann = [{"id": "axes", "kind": "triad", "at": [-250, -250, 0], "length": 330}]

# the sinusoidal transverse load: n x n arrows on the plate's upper face, the
# length of each arrow growing with q at its foot, so the arrow field is the load
# surface itself; l0 is a common floor so that even the corner arrows of a
# vanishing q keep a visible tail
for i in range(n):
    x = (2*i + 1)*L/(2*n)
    for j in range(n):
        y = (2*j + 1)*L/(2*n)
        q = math.sin(math.pi*x/L)*math.sin(math.pi*y/L)
        ann.append({"id": f"q_{i}{j}", "kind": "force", "at": [x, y, t], "dir": [0, 0, -1],
                    "length": l0 + l1*q, "head": True})

ann += [
    {"id": "load", "kind": "label", "at": [L, L, t], "role": "load",
     "text": r"$\bar q=\sin(\pi x)\sin(\pi y)$", "screen_offset": [0, 215]},
    {"id": "Lx", "kind": "dim", "p0": [0, 0, t/2], "p1": [L, 0, t/2],
     "offset": [0, -1, 0], "dist": 170, "text": "$L$"},
    {"id": "Ly", "kind": "dim", "p0": [0, 0, t/2], "p1": [0, L, t/2],
     "offset": [-1, 0, 0], "dist": 170, "text": "$L$"},
]


def model_1d(tag):
    """The 3D point (along-axis position a, height z) and the axis of a 1D model."""
    if tag == "x":
        return (lambda a, z: [a, L + d1, z]), [1, 0, 0]
    return (lambda a, z: [L + d1, a, z]), [0, 1, 0]


def seg(prefix, p0, p1, color="0.1", lw=1.2):
    """One straight piece of the line work, drawn on the scene's camera."""
    ann.append({"id": prefix, "kind": "line", "style": "solid", "color": color,
                "lw": lw, "p0": list(p0), "p1": list(p1)})


def ring(prefix, c, r, u, n=14, color="0.1", lw=0.9):
    """A small circle in the vertical plane through u: a hinge or a roller."""
    u = np.asarray(u, float) / np.linalg.norm(u)
    c = np.asarray(c, float)
    pts = [c + r*(math.cos(a)*u + math.sin(a)*np.array([0.0, 0.0, 1.0]))
           for a in (2*math.pi*k/n for k in range(n + 1))]
    for k in range(n):
        seg(f"{prefix}{k}", pts[k], pts[k + 1], color, lw)


for tag in ("x", "y"):
    pt, u = model_1d(tag)
    # the beam line on the mid-surface, its arrows onto the line, and the sine
    # envelope through their tails
    seg(f"beaml_{tag}", pt(r_h + 12, t/2), pt(L - r_h - 12, t/2), "0.1", 1.6)
    for i in range(m):
        a = (2*i + 1)*L/(2*m)
        ann.append({"id": f"q1d_{tag}{i}", "kind": "force", "at": pt(a, t/2),
                    "dir": [0, 0, -1], "head": True,
                    "length": b0 + b1*math.sin(math.pi*a/L)})
    curve = [(a, t/2 + b0 + b1*math.sin(math.pi*a/L)) for a in
             (k*L/20 for k in range(21))]
    for k in range(len(curve) - 1):
        seg(f"env_{tag}{k}", pt(*curve[k]), pt(*curve[k + 1]), sk.LOAD, 1.3)
    ann.append({"id": f"lab_{tag}", "kind": "label", "role": "load",
                "at": pt(L/2, t/2 + b0 + b1 + 150),
                "text": rf"$\bar q=\sin(\pi {tag})$"})

    # the two supports: a white hinge disc at each beam end, a triangle down to
    # the filled ground bar at the pinned end, and the same triangle on two
    # rollers at the other end
    nrm = np.cross(np.asarray(u, float), [0, 0, 1])      # across the drawing plane
    for end, z_base in ((0, z_base_pin), (L, z_base_rol)):
        c = np.asarray(pt(end, t/2), float)
        parts.append({"id": f"hinge_{tag}{end}", "kind": "rod", "color": "white",
                      "args": {"p0": list(c - dis_t/2*nrm), "p1": list(c + dis_t/2*nrm),
                               "r": r_h}})
        apex, left, right = (pt(end, z_apex), pt(end - w_tri/2, z_base),
                             pt(end + w_tri/2, z_base))
        for k, (p, q) in enumerate(((left, right), (right, apex), (apex, left))):
            seg(f"tri_{tag}{end}_{k}", p, q, "0.1", 1.2)
        # the ground, filled as a bar rather than hatched
        box_along(f"ground_{tag}{end}", pt(end, 0.0), u, nrm, 2*0.95*w_tri, gnd_t,
                  g1 - gnd_h, g1, "block")
    for k in (-1, 1):
        ring(f"roller_{tag}{k}", pt(L + k*w_tri/4, g1 + r_r), r_r, u)

spec = {
    "name": "sinusoidal_plate",
    "params": {"L": L, "t": t, "r_h": r_h, "w_tri": w_tri, "n": n, "l0": l0,
               "l1": l1, "d1": d1},
    "parts": parts,
    "annotations": ann,
    "cameras": [{"az": -35, "el": 32}, {"az": -40, "el": 30}, {"az": -45, "el": 30}],
    "pick": {"az": -35, "el": 32},
    "light": [-0.25, 0.45, 0.85],
    "width_in": 5.2,
}

# ------------------------------------------------------------------ render
OUT.mkdir(parents=True, exist_ok=True)
(OUT/"sinusoidal_plate.json").write_text(json.dumps(spec, indent=1))
report = s3.run(spec, OUT)
print("physics:", report["physics"])
for c in report["candidates"]:
    print(f"  az {c['camera']['az']:>4} el {c['camera']['el']:>3}  score {c['score']:.2f}"
          f"  occluded {sum(c['occluded_fraction'].values()):.2f}"
          f"  overlaps {len(c['label_overlaps'])}")
print("used:", report["used"])

# The figure of the paper: the PNG preview (200 dpi) and a PDF for XeLaTeX.
# The PDF is the SVG converted by rsvg-convert, so it stays vector; the text is
# already paths (matplotlib's svg.fonttype), so no font is needed to view it.
# Labels that sit outside the frame of the solids are listed by the report under
# labels_outside; bbox_inches='tight' keeps them in the saved figure.
shutil.copy(OUT/"sinusoidal_plate.png", FIGDIR/"sinusoidal_plate.png")
subprocess.run(["rsvg-convert", "-f", "pdf", "-o", str(FIGDIR/"sinusoidal_plate.pdf"),
                str(OUT/"sinusoidal_plate.svg")], check=True)
