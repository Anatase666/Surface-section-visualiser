# Plane–Surface Cross-Section Visualizer

An interactive 3D tool for exploring how a **cutting plane** intersects a
**3D surface** `z = f(x, y)`. A [VPython](https://vpython.org/) window shows
the surface and a movable/tiltable plane in 3D, while a linked
[Matplotlib](https://matplotlib.org/) window plots the *true* cross-section
curve in real time as you drag three sliders (tilt, X shift, Y shift).

The cross-section is computed as an actual **implicit intersection curve**
(`z_plane(u, v) − f(x(u, v), y(u, v)) = 0`, extracted with a contour
algorithm), not a naive "height profile along a line" — so closed curves
(ellipses, circles) and open branches all come out correctly, matching what
you see rendered in the 3D view.

## Why this exists / applications

This project was built as a teaching/exploration aid for anyone who needs an
intuitive feel for "what does a slice through a 3D surface actually look
like", including:

- **Geology** — building intuition for how geological **cross-sections**
  (a vertical "cut" through folded/tilted rock strata) relate to the 3D
  shape of the layers/terrain they're cut from. Load a real or synthetic
  elevation/horizon surface and practice reading how the section shape
  changes as the cutting plane's position and dip (tilt) change.
- **Mathematics** — visualizing intersections of a plane with a surface
  (e.g. conic-section-style curves cut from a paraboloid, saddle, or any
  custom function), as a hands-on companion to analytic geometry /
  multivariable calculus.
- **Engineering / technical (descriptive) graphics** — constructing the
  *true section* of a solid or surface cut by an oblique plane, a classic
  descriptive-geometry exercise, here made interactive instead of manual.

## Features

- Two ways to define the 3D surface:
  - **Function mode** — any symbolic `z = f(x, y)` expression (via SymPy).
  - **File mode** — load a real point cloud (`x, y, z` per line, CSV or
    whitespace-separated) and the surface is reconstructed via
    interpolation. Points do **not** need to lie on a regular grid.
- Three sliders control the cutting plane in real time:
  - **Plane angle** — tilts the plane about its hinge line.
  - **Shift X / Shift Y** — translates the plane across the surface.
- The 3D plane patch is drawn intentionally larger than the surface's own
  footprint, so it's always visible extending past the surface's edges.
- The cross-section is recomputed as a true implicit curve on every slider
  move and plotted in the plane's own undistorted `(u, v)` coordinates
  (equal aspect ratio), so shapes like circles/ellipses look correct rather
  than stretched.

## Requirements & installation

Requires **Python 3.9+**.

```bash
# 1. clone the repository
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>

# 2. create and activate a virtual environment
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate

# 3. install dependencies
pip install -r requirements.txt
```

An equivalent `environment.yml` is provided for Conda users:

```bash
conda env create -f environment.yml
conda activate plane-section-viz
```

> **Note:** VPython renders its 3D scene in a browser tab/window (it uses a
> local WebSocket server under the hood), so running the script will open
> your default browser automatically — this is expected, not an error.

## Usage

```bash
python main.py
```

Two windows/tabs will open:

1. A **browser tab** with the 3D scene (surface + cutting plane).
2. A **Matplotlib window** with three sliders and the live cross-section
   plot.

### Switching between function mode and file mode

Near the top of `main.py`:

```python
SURFACE_MODE = "file"          # "function"  or  "file"
FILE_PATH = "example_surface_points.csv"   # used only when SURFACE_MODE == "file"
```

- Set `SURFACE_MODE = "function"` and edit the `function_expr = x**2 + y**2`
  line to visualize any analytic surface instead.
- Set `SURFACE_MODE = "file"` and point `FILE_PATH` at your own point-cloud
  file to visualize real/custom surface data. The working domain (plane
  bounds, slider ranges) is automatically taken from the extent of the
  loaded points.

### Point-cloud file format

A plain text (`.csv` or `.txt`) file, **one point per line**:

```
x,y,z
-2.500,-2.500,-0.0284
-2.500,-2.300,-0.0585
...
```

- Values may be separated by **commas or whitespace** (auto-detected).
- An optional header row (e.g. `x,y,z`) is fine — non-numeric lines are
  skipped automatically.
- Points do not need to be on a regular grid; a scattered cloud works too.

A ready-to-use example, `example_surface_points.csv` (676 points, a
synthetic undulating surface reminiscent of folded rock strata), is
included in this repository so you can try file mode immediately.

## Project structure

```
.
├── main.py                       # the application
├── example_surface_points.csv    # sample point-cloud surface (file mode demo)
├── requirements.txt              # pip dependencies
├── environment.yml               # optional Conda environment spec
├── .gitignore
├── LICENSE
└── README.md
```

## Known limitations

- Very large point clouds (tens of thousands of points) will make the
  scattered-data interpolation and the vpython triangle-mesh construction
  noticeably slower; a few hundred to a few thousand points works smoothly.
- The cross-section contour is computed on a fixed-resolution grid
  (`n_grid` in `get_section_contours`); increase it for a smoother curve at
  the cost of speed, or decrease it if sliders feel laggy.

## License

Released under the [MIT License](LICENSE) — see the file for details.
