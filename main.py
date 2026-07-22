from vpython import *
from math import *
import sympy as sp
import numpy as np
import csv

from scipy.interpolate import griddata, LinearNDInterpolator

import matplotlib
matplotlib.use("TkAgg")
# Force an interactive backend explicitly, and do it AFTER importing vpython:
# some environments (e.g. PyCharm's built-in "SciView" tool window, a plain
# 'Agg' backend, or vpython itself resetting the backend on import) render
# static images only and silently ignore widget events, so the Slider looks
# "dead" even though the callback code itself is correct. TkAgg ships with
# the standard Python installation, so it works without extra installs;
# swap for "QtAgg" if you already have PyQt/PySide installed and prefer that.
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

print("Matplotlib backend in use:", matplotlib.get_backend())


def cos_angle_btw_vec(vec1, vec2):
    return dot(vec1, vec2)/(vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z)


def newthsolve(equation, x0=0, dx=0.01, max_iter=100):

    x = sp.Symbol('x')
    equation = equation + x
    derivative = equation.diff(x)
    iter = 0
    while iter < max_iter:
        iter += 1
        expr = x - equation/derivative
        x0 = expr.evalf(subs={x: x0})
        if abs(x0) < dx:
            return x0, iter
    return None, iter


def newtfursolve(equation, a=0, b=0, dx=0.01, max_iter=100):

    x = sp.Symbol('x')
    equation = equation + x
    derivative = equation.diff(x)
    iter = 0
    while iter < max_iter:
        iter += 1
        expr = x - equation/derivative
        a = expr.evalf(subs={x: a})
        b = expr.evalf(subs={x: b})
        if abs(a-b) < dx:
            return a, iter
    return None, iter


class Line:
    def __init__(self, direct=None, spot=None):
        if direct is None:
            direct = vector(0, 0, 0)
        if spot is None:
            spot = vector(0, 0, 0)
        self.direct = norm(direct)
        self.spot = spot


class Plane:
    def __init__(self, line=None, angle=None, normal=None, spot=None):

        if normal is not None and spot is not None:
            self.spot = spot
            self.normal = normal
            self.angle = diff_angle(normal, vector(0, 0, 1))
            self.line = Line()
            self.line.spot = spot
            self.line.direct = norm(cross(normal, vector(0, 0, 1)))

        elif line is not None and angle is not None:
            self.line = line
            self.spot = vector(0, 0, 0)
            self.spot = line.spot
            self.normal = norm(rotate(vector(0, 0, 1), angle=angle, axis=line.direct))
        else:
            raise TypeError('Plane is not initialised')

        self.direct_vec1 = rotate(self.normal, angle=pi/2, axis=self.line.direct)
        self.direct_vec2 = norm(cross(vector(0, 0, 1), self.normal))
        # two complanar vectors
        self.d = -(self.normal.x*self.spot.x + self.normal.y*self.spot.y + self.normal.z*self.spot.z)

        # triangles created by draw() are kept here so a previous drawing of
        # this (or a replaced) plane can be hidden/removed before redrawing
        self.tris = []

    def calc_coord_from_param(self, p1, p2):
        # calculate coordinates of plane's spots using two params
        x1 = self.spot.x + p1 * self.direct_vec1.x + p2 * self.direct_vec2.x
        y1 = self.spot.y + p1 * self.direct_vec1.y + p2 * self.direct_vec2.y
        z1 = self.spot.z + p1 * self.direct_vec1.z + p2 * self.direct_vec2.z
        return x1, y1, z1

    def calc_borders_for_params(self, xmin, ymin, xmax, ymax):
        b_a = self.direct_vec2.y/self.direct_vec2.x
        m_min = (ymin - self.spot.y - b_a*(xmin - self.spot.x))/(self.direct_vec1.y - b_a*self.direct_vec1.x)
        n_min = (xmin - self.spot.x - self.direct_vec1.x*m_min)/self.direct_vec2.x

        m_max = (ymax - self.spot.y - b_a * (xmax - self.spot.x)) / (self.direct_vec1.y - b_a * self.direct_vec1.x)
        n_max = (xmax - self.spot.x - self.direct_vec1.x * m_max) / self.direct_vec2.x
        return m_min, n_min, m_max, n_max

    def calc_equation_of_plane(self):
        # gets symbolic system of equations for plane
        m, n = sp.symbols('m n')
        x1 = self.spot.x + m * self.direct_vec1.x + n * self.direct_vec2.x
        y1 = self.spot.y + m * self.direct_vec1.y + n * self.direct_vec2.y
        z1 = self.spot.z + m * self.direct_vec1.z + n * self.direct_vec2.z
        return x1, y1, z1

    def draw(self, xmin, ymin, xmax, ymax, extent_scale=1.0):
        # Build an actual rectangle (4 distinct corners). extent_scale lets
        # the plane patch be drawn larger than the surface's own bounding
        # box (e.g. 1.4x) so the cutting plane is always clearly visible
        # sticking out past the edges of the surface, instead of stopping
        # exactly at them.
        half_size = sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2) / 2 * extent_scale

        c1 = self.spot + half_size * self.direct_vec1 + half_size * self.direct_vec2
        c2 = self.spot + half_size * self.direct_vec1 - half_size * self.direct_vec2
        c3 = self.spot - half_size * self.direct_vec1 - half_size * self.direct_vec2
        c4 = self.spot - half_size * self.direct_vec1 + half_size * self.direct_vec2

        vert1 = vertex(pos=c1, color=vector(1, 1, 1))
        vert2 = vertex(pos=c2, color=vector(1, 1, 1))
        vert3 = vertex(pos=c3, color=vector(1, 1, 1))
        vert4 = vertex(pos=c4, color=vector(1, 1, 1))

        self.tris = [
            triangle(vs=[vert1, vert2, vert3]),
            triangle(vs=[vert1, vert3, vert4]),
        ]
        return self.tris

    def clear(self):
        # hides the triangles this plane previously drew in the 3D scene
        for t in self.tris:
            t.visible = False
        self.tris = []


def dist_from_spot_to_plane(spot, plane):
    return abs((spot.pos.x*plane.normal.x + spot.pos.y*plane.normal.y + spot.pos.z*plane.normal.z + plane.d))\
           / sqrt(plane.normal.x**2 + plane.normal.y**2 + plane.normal.z**2)


def rangein(x1, x2, dx):
    # makes a list of values with equal gaps
    L = []
    L.append(x1)
    xv = x1
    while xv < x2:
        xv += dx
        L.append(xv)
    return L


def FindIn2D(LV, dot, dx):
    # finds discrete coordinates of spot in 2D list
    length = len(LV)
    for nx in range(length):
        for ny in range(length):
            if abs(nx.pos.x-dot[0]) <= dx and abs(ny.pos.y-dot[1]) <= dx:
                return nx, ny


# ---------------------------------------------------------------------------
# Surface sources -- "function" (analytic z = f(x, y)) or "file" (a cloud of
# (x, y, z) points loaded from disk). Both end up as the same thing: a fast
# vectorised numpy callable f_num(X, Y) -> Z, which is all the rest of the
# program (3D grid builder + matplotlib section) actually needs.
# ---------------------------------------------------------------------------

def make_function_surface(expr, x_sym, y_sym):
    """Turns a sympy expression into a vectorised numpy callable z = f(x, y)."""
    return sp.lambdify((x_sym, y_sym), expr, 'numpy')


def load_points_file(path):
    """
    Reads a text/CSV file describing a surface as a cloud of 3D points.

    Format: one point per line, three numbers (x, y, z), separated by
    commas or whitespace. An optional header line (e.g. "x,y,z") is
    detected and skipped automatically. Points do NOT need to lie on a
    regular grid -- a Delaunay-based linear interpolator is built from
    them, so the surface can be resampled onto whatever (x, y) grid the
    rest of the program needs (this also means a coarse or irregular
    point cloud is fine; it doesn't have to match `dx` below).

    Returns: f_num, xmin, xmax, ymin, ymax  (domain = extent of the data)
    """
    pts = []
    with open(path, 'r', newline='') as f:
        sample = f.read(2048)
        f.seek(0)
        delimiter = ',' if ',' in sample else None
        rows = csv.reader(f, delimiter=delimiter) if delimiter else (line.split() for line in f)
        for row in rows:
            if not row:
                continue
            try:
                x1, y1, z1 = float(row[0]), float(row[1]), float(row[2])
            except (ValueError, IndexError):
                continue  # skips header row / blank / malformed lines
            pts.append((x1, y1, z1))

    if len(pts) < 3:
        raise ValueError(f"Not enough valid (x, y, z) points found in {path}")

    pts = np.array(pts)
    xs, ys, zs = pts[:, 0], pts[:, 1], pts[:, 2]
    interpolator = LinearNDInterpolator(list(zip(xs, ys)), zs)

    def f_num(X, Y):
        Z = interpolator(X, Y)
        Z = np.asarray(Z, dtype=float)
        # points outside the convex hull of the loaded cloud come back as
        # NaN from the linear interpolator -- fall back to nearest-neighbour
        # so the surface/section stay well defined all the way to the edges
        if np.any(np.isnan(Z)):
            nn = griddata((xs, ys), zs, (X, Y), method='nearest')
            Z = np.where(np.isnan(Z), nn, Z)
        return Z

    return f_num, float(xs.min()), float(xs.max()), float(ys.min()), float(ys.max())


def eval_z_scalar(f_num, x1, y1):
    """Evaluates the vectorised surface function at a single (x, y) point."""
    z = f_num(np.array([x1]), np.array([y1]))
    return float(np.asarray(z).reshape(-1)[0])


def model_func(f_num, xmin, ymin, xmax, ymax, dx):
    """
    Builds the vpython grid of vertices/triangles for a surface z = f(x, y).

    Note: colouring and triangle generation now happen ONCE, after the full
    (x, y) grid has been sampled -- in the original version they were
    re-run from scratch on every step of the outer x-loop (an accidental
    O(n^3) redraw). That was barely noticeable on the old 2x2 domain, but
    with the enlarged domain below it's the difference between an instant
    redraw and a very long wait, so it's fixed here.
    """
    xs_range = rangein(xmin, xmax, dx)
    ys_range = rangein(ymin, ymax, dx)

    LV = []
    zmax = -1e100
    zmin = 1e100
    for x1 in xs_range:
        TempLV = []
        for y1 in ys_range:
            z1 = eval_z_scalar(f_num, x1, y1)
            zmax = max(zmax, z1)
            zmin = min(zmin, z1)
            TempLV.append(vertex(pos=vector(x1, y1, z1), color=vector(1, 0.1, 1)))
        LV.append(TempLV)

    N = len(LV)
    M = len(LV[0]) if N else 0
    zrange = zmax if zmax != 0 else 1.0

    for ix in range(N):
        for iy in range(M):
            blue = LV[ix][iy].pos.z / zrange
            if blue > 0:
                green = 1 - blue
                red = 0
            else:
                green = blue
                blue = 1 - blue
                red = 0
            LV[ix][iy].color = vector(red, green, blue)

    for ix in range(1, N - 1):
        for iy in range(1, M - 1):
            triangle(vs=[LV[ix][iy], LV[ix + 1][iy], LV[ix][iy + 1]])
            triangle(vs=[LV[ix][iy], LV[ix - 1][iy], LV[ix][iy + 1]])
            triangle(vs=[LV[ix][iy], LV[ix + 1][iy], LV[ix][iy - 1]])
            triangle(vs=[LV[ix][iy], LV[ix - 1][iy], LV[ix][iy - 1]])

    return zmin, zmax, LV


# ---------------------------------------------------------------------------
# Parametric section computation (unchanged logic, now takes f_num directly
# instead of a sympy expression + symbols, so it works for both surface
# sources without caring which one is active).
# ---------------------------------------------------------------------------

_contour_fig = Figure()
FigureCanvasAgg(_contour_fig)
_contour_ax = _contour_fig.add_subplot(111)


def get_section_contours(plane, f_num, half_extent, n_grid=400):
    u = np.linspace(-half_extent, half_extent, n_grid)
    v = np.linspace(-half_extent, half_extent, n_grid)
    U, V = np.meshgrid(u, v)

    Xs = plane.spot.x + U * plane.direct_vec1.x + V * plane.direct_vec2.x
    Ys = plane.spot.y + U * plane.direct_vec1.y + V * plane.direct_vec2.y
    Zplane = plane.spot.z + U * plane.direct_vec1.z + V * plane.direct_vec2.z

    G = Zplane - f_num(Xs, Ys)

    _contour_ax.clear()
    cs = _contour_ax.contour(U, V, G, levels=[0.0])
    segments = [seg.copy() for seg in cs.allsegs[0]]
    return segments


# ---------------------------------------------------------------------------
# Surface configuration -- pick ONE of the two modes below.
# ---------------------------------------------------------------------------
print("Choose a mode: \n (1) - file \n (2) - function")
SURFACE_MODE = input()         # "function" -- analytic z = f(x, y) (sympy)
                                    # "file"     -- load (x, y, z) points from disk

if SURFACE_MODE == '1':
    print("Please, type file's path")
    FILE_PATH = "example_surface_points.csv"   # name default
    f = input()
    if f:
        FILE_PATH = f

# Domain of the surface/plane. Enlarged from the original -1..1 square so
# there is visibly more room for the cutting plane to move/tilt in.
xmin, ymin, xmax, ymax = -2.5, -2.5, 2.5, 2.5
dx = 0.1                           # sampling step for the vpython surface
PLANE_EXTENT_SCALE = 1.4           # plane patch drawn 1.4x larger than the
                                    # surface's bounding box, so it visibly
                                    # extends past the surface's edges

x, y = sp.symbols("x y")

if SURFACE_MODE == '2':
    print("Type function")
    func = input()
    if func:
        x, y = sp.symbols('x y')
        function_expr = sp.sympify(func)
    else:
        function_expr = x**2 + y**2       # function default

    f_num = make_function_surface(function_expr, x, y)

elif SURFACE_MODE == 1:
    f_num, fxmin, fxmax, fymin, fymax = load_points_file(FILE_PATH)
    # the working domain becomes the extent of the loaded point cloud
    xmin, xmax, ymin, ymax = fxmin, fxmax, fymin, fymax

else:
    raise ValueError("SURFACE_MODE must be '1' or '2'")

# 3D surface in vpython
zmin, zmax, LV = model_func(f_num, xmin, ymin, xmax, ymax, dx)

# base spot/direction of the cutting line: spot.x and spot.y are the
# quantities the two shift sliders move; direct stays fixed so that only
# position (x, y sliders) and tilt (angle slider) change the plane.
BASE_Z = 0.3
line1 = Line(vector(1, 1, 0), vector(0, 0, BASE_Z))
plane1 = Plane(normal=vector(1, 1, 0), spot=vector(0, 0, 0))
plane1.draw(xmin, ymin, xmax, ymax, extent_scale=PLANE_EXTENT_SCALE)

# half-size of the plane patch actually drawn in 3D (same formula as
# Plane.draw uses, including the extent scale) -- using the same extent for
# the (u, v) contour grid means the matplotlib section is computed over
# exactly the white patch you see in the vpython window, so the two views
# always agree.
half_extent = sqrt((xmax - xmin) ** 2 + (ymax - ymin) ** 2) / 2 * PLANE_EXTENT_SCALE

segments0 = get_section_contours(plane1, f_num, half_extent)

# ---------------------------------------------------------------------------
# Interactive matplotlib plot: three sliders control the cutting plane --
# "Plane angle" rotates it around line1, "Shift X"/"Shift Y" translate it.
# The plot shows the section in the plane's OWN 2D coordinates (u along
# direct_vec1, v along direct_vec2); since these two basis vectors are
# orthonormal, an equal aspect ratio keeps the true shape undistorted.
# ---------------------------------------------------------------------------

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.32)

ax.set_xlabel("u (in-plane, along direct_vec1)")
ax.set_ylabel("v (in-plane, along direct_vec2)")
ax.set_aspect("equal", adjustable="box")
ax.grid(True)

section_artists = []


def draw_segments(segments):
    global section_artists
    for artist in section_artists:
        artist.remove()
    section_artists = []
    for seg in segments:
        (line_artist,) = ax.plot(seg[:, 0], seg[:, 1], "-", color="blue", label="Section")
        section_artists.append(line_artist)
    ax.set_xlim(-half_extent, half_extent)
    ax.set_ylim(-half_extent, half_extent)
    if segments:
        ax.legend(["Section"])
    else:
        ax.legend([])


draw_segments(segments0)

ax_angle = plt.axes([0.2, 0.17, 0.6, 0.03])
angle_slider = Slider(ax_angle, "Plane angle", 0, pi, valinit=0)

ax_shift_x = plt.axes([0.2, 0.10, 0.6, 0.03])
shift_x_slider = Slider(ax_shift_x, "Shift X", xmin, xmax, valinit=0)

ax_shift_y = plt.axes([0.2, 0.03, 0.6, 0.03])
shift_y_slider = Slider(ax_shift_y, "Shift Y", ymin, ymax, valinit=0)


def update(val):
    global plane1
    try:
        # move the cutting line to the (x, y) chosen by the two shift
        # sliders, keeping the same fixed direction and height as before
        line1.spot = vector(shift_x_slider.val, shift_y_slider.val, BASE_Z)

        new_plane = Plane(line=line1, angle=angle_slider.val)

        # remove the previously drawn plane patch and draw the new one, so
        # the 3D tilt/position actually updates together with the slider(s)
        plane1.clear()
        new_plane.draw(xmin, ymin, xmax, ymax, extent_scale=PLANE_EXTENT_SCALE)
        plane1 = new_plane

        segments = get_section_contours(plane1, f_num, half_extent)
        draw_segments(segments)
        fig.canvas.draw_idle()
        fig.canvas.flush_events()
        print("angle:", angle_slider.val, "shift_x:", shift_x_slider.val, "shift_y:", shift_y_slider.val)
    except Exception:
        import traceback
        with open("main_slider_error.log", "w") as f:
            traceback.print_exc(file=f)
        print("ERROR - see main_slider_error.log")


angle_slider.on_changed(update)
shift_x_slider.on_changed(update)
shift_y_slider.on_changed(update)

plt.show()