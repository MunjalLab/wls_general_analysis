"""
Geodesic-band centerline extraction for a hollow, open tube-wall mesh in Blender.

What this script is for
-----------------------
This script is designed for a 3D segmentation mesh that is shaped like a hollow,
curved tube. The mesh should be the surface/wall of the tube, not a filled solid.
The tube should have two open ends and no caps.

Before running, check the mesh in Blender:
    Edit Mode > Select > Select All by Trait > Non-Manifold

Ideally, the only selected non-manifold edges should be the two open rims. Fill
any accidental small holes in the tube wall before using this script.

Method overview
---------------
This script avoids unstable plane-slicing. Instead, it uses a geodesic coordinate
along the tube surface:

    1. Detect the two open boundary rim loops.
    2. Build a weighted graph from mesh edges.
    3. Compute shortest-path surface distance from every vertex to rim 1.
    4. Compute shortest-path surface distance from every vertex to rim 2.
    5. Define a normalized longitudinal coordinate:

           s = distance_to_rim_1 / (distance_to_rim_1 + distance_to_rim_2)

       s is near 0 at rim 1, near 1 at rim 2, and near 0.5 around the middle.
    6. Group vertices into thin geodesic bands by similar s values.
    7. For each band, estimate a local cross-section center and radius.
    8. Smooth the sequence of centers to produce the final centerline.

Outputs created in Blender
--------------------------
    geodesic_band_midline          A curve object showing the centerline.
    geodesic_band_cross_sections   Optional visual radius circles.
    tube_end_1_centroid            Marker at one open rim centroid.
    tube_end_2_centroid            Marker at the other open rim centroid.

CSV outputs
-----------
If SAVE_SECTION_CSV is True, the script writes one row per centerline point:
    geodesic_band_centerline_sections.csv

If SAVE_SUMMARY_CSV is True, the script writes the console-style summary metrics:
    geodesic_band_summary_metrics.csv

If the .blend file is saved, CSV files are saved next to the .blend file.
If the .blend file is not saved, CSV files are saved to your Desktop.

How to run
----------
    1. Select the tube mesh object in Object Mode.
    2. Open this script in Blender's Scripting workspace.
    3. Press Run Script.

The results are also stored in the Python variable named `results`.
"""

import bpy
import numpy as np
import math
import heapq
import csv
import os


# ============================================================
# USER SETTINGS
# ============================================================

# Number of points sampled along the final centerline.
# Increase for more detailed output. Decrease if your mesh is sparse.
N_SECTIONS = 100

# Smoothing applied to the final sequence of center points.
# Increase these if the centerline is slightly jagged.
# Decrease these if the centerline cuts corners too much.
SMOOTHING_PASSES = 8
SMOOTHING_STRENGTH = 0.55

# Blender visualization output toggles.
CREATE_CENTERLINE_OBJECT = True
CREATE_SECTION_CIRCLES = True
CREATE_ENDPOINT_MARKERS = True

# CSV output toggles.
SAVE_SECTION_CSV = True
SAVE_SUMMARY_CSV = True

# Object and file names.
CENTERLINE_OBJECT_NAME = "geodesic_band_midline"
SECTION_COLLECTION_NAME = "geodesic_band_cross_sections"
SECTION_CSV_FILENAME = "geodesic_band_centerline_sections.csv"
SUMMARY_CSV_FILENAME = "geodesic_band_summary_metrics.csv"

# Visual thickness of generated curve objects.
# Adjust these depending on the scale of your mesh.
CENTERLINE_BEVEL_DEPTH = 0.01
SECTION_CIRCLE_BEVEL_DEPTH = 0.004
ENDPOINT_MARKER_RADIUS = 0.05

# Minimum number of mesh vertices required in each longitudinal band.
# If your mesh is sparse and some bands fail, reduce this to 8 or so.
MIN_POINTS_PER_SECTION = 12

# If too few vertices are found for a band, the script expands the band width.
MAX_BIN_EXPANSIONS = 8
BIN_EXPANSION_FACTOR = 1.5

# If True, radius is computed as the mean distance from band vertices to the
# fitted section center. This is usually stable for biological segmentation meshes.
# If False, radius is taken directly from the 2D least-squares circle fit.
RADIUS_AS_MEAN_DISTANCE = True

# Optional radius sanity limits. Leave as None unless you know the expected range.
MIN_REASONABLE_RADIUS = None
MAX_REASONABLE_RADIUS = None


# ============================================================
# GENERAL HELPERS
# ============================================================

def normalize(v, eps=1e-12):
    """Return a unit-length version of vector v."""
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n < eps:
        return v * 0.0
    return v / n


def get_selected_mesh_object():
    """Get the currently selected Blender object and make sure it is a mesh."""
    obj = bpy.context.object
    if obj is None or obj.type != "MESH":
        raise RuntimeError("Please select a mesh object before running the script.")
    return obj


def get_world_vertices_faces_edges(obj):
    """
    Extract evaluated mesh vertices, faces, and unique edges in world coordinates.

    This uses the evaluated object, so active modifiers may affect the topology.
    If boundary-loop counts are unexpected, try applying or disabling modifiers.
    """
    depsgraph = bpy.context.evaluated_depsgraph_get()
    eval_obj = obj.evaluated_get(depsgraph)
    mesh = eval_obj.to_mesh()

    vertices = np.array(
        [eval_obj.matrix_world @ v.co for v in mesh.vertices],
        dtype=float
    )

    faces = []
    edges_set = set()

    for poly in mesh.polygons:
        face = list(poly.vertices)
        faces.append(face)

        for i in range(len(face)):
            a = face[i]
            b = face[(i + 1) % len(face)]
            edges_set.add(tuple(sorted((a, b))))

    edges = np.array(sorted(edges_set), dtype=int)

    eval_obj.to_mesh_clear()
    return vertices, faces, edges


def polyline_length(points):
    """Compute total length of a 3D polyline."""
    points = np.asarray(points, dtype=float)
    if len(points) < 2:
        return 0.0
    segment_lengths = np.linalg.norm(points[1:] - points[:-1], axis=1)
    return float(np.sum(segment_lengths))


def get_output_directory():
    """
    Return a writable output directory.

    If the .blend file is saved, output files are saved next to it.
    If not, output files are saved to the Desktop.
    """
    if bpy.data.filepath:
        return os.path.dirname(bpy.data.filepath)
    return os.path.join(os.path.expanduser("~"), "Desktop")


# ============================================================
# BOUNDARY LOOP DETECTION
# ============================================================

def find_boundary_loops(faces):
    """
    Find connected components of boundary edges.

    For a clean hollow open tube wall, this should return exactly two loops:
        loop 1 = one open rim
        loop 2 = other open rim
    """
    edge_count = {}

    # Count how many faces use each edge.
    # Boundary edges are edges used by exactly one face.
    for face in faces:
        for i in range(len(face)):
            a = face[i]
            b = face[(i + 1) % len(face)]
            edge = tuple(sorted((a, b)))
            edge_count[edge] = edge_count.get(edge, 0) + 1

    boundary_edges = [edge for edge, count in edge_count.items() if count == 1]

    if len(boundary_edges) == 0:
        raise RuntimeError("No boundary edges found. Expected two open tube rims.")

    # Build adjacency graph containing only boundary edges.
    adjacency = {}
    for a, b in boundary_edges:
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)

    # Connected components of this graph are boundary loops.
    visited = set()
    loops = []

    for start in adjacency.keys():
        if start in visited:
            continue

        stack = [start]
        component = []
        visited.add(start)

        while stack:
            v = stack.pop()
            component.append(v)

            for nb in adjacency[v]:
                if nb not in visited:
                    visited.add(nb)
                    stack.append(nb)

        loops.append(component)

    loops = sorted(loops, key=len, reverse=True)
    return loops


# ============================================================
# GEODESIC DISTANCE ON THE MESH SURFACE
# ============================================================

def build_weighted_vertex_graph(vertices, edges):
    """
    Build a weighted graph from mesh edges.

    Each mesh vertex is a graph node, and each mesh edge is weighted by its
    3D length. Dijkstra's algorithm on this graph approximates geodesic distance
    along the mesh surface.
    """
    graph = [[] for _ in range(len(vertices))]

    for a, b in edges:
        weight = float(np.linalg.norm(vertices[a] - vertices[b]))
        graph[a].append((b, weight))
        graph[b].append((a, weight))

    return graph


def dijkstra_from_sources(graph, sources):
    """
    Compute shortest-path distance from a set of source vertices to all vertices.

    In this script, the sources are all vertices in one open rim loop.
    """
    n = len(graph)
    dist = np.full(n, np.inf, dtype=float)
    heap = []

    for source in sources:
        source = int(source)
        dist[source] = 0.0
        heapq.heappush(heap, (0.0, source))

    while heap:
        d, v = heapq.heappop(heap)

        # Skip stale heap entries.
        if d != dist[v]:
            continue

        for nb, weight in graph[v]:
            new_dist = d + weight
            if new_dist < dist[nb]:
                dist[nb] = new_dist
                heapq.heappush(heap, (new_dist, nb))

    return dist


# ============================================================
# LOCAL CENTER AND RADIUS ESTIMATION
# ============================================================

def fit_pca_basis(points):
    """
    Fit a local PCA basis to a point cloud.

    Returns centroid and three orthogonal principal axes. For each geodesic band,
    the points should form a ring-like section. PCA gives a stable local frame
    for fitting a circle.
    """
    points = np.asarray(points, dtype=float)
    centroid = points.mean(axis=0)
    centered = points - centroid
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    return centroid, vh[0], vh[1], vh[2]


def fit_circle_2d(xy):
    """
    Least-squares circle fit in 2D.

    Fits:
        x^2 + y^2 + A*x + B*y + C = 0

    Returns the 2D circle center and radius.
    """
    xy = np.asarray(xy, dtype=float)
    if len(xy) < 3:
        raise RuntimeError("Need at least 3 points to fit a circle.")

    x = xy[:, 0]
    y = xy[:, 1]

    M = np.column_stack([x, y, np.ones_like(x)])
    rhs = -(x * x + y * y)

    A, B, C = np.linalg.lstsq(M, rhs, rcond=None)[0]

    center = np.array([-A / 2.0, -B / 2.0])
    radius_squared = center[0] * center[0] + center[1] * center[1] - C

    if radius_squared <= 0:
        raise RuntimeError("Circle fit failed with non-positive radius.")

    radius = math.sqrt(radius_squared)
    return center, radius


def robust_section_center_and_radius(points):
    """
    Estimate the center and radius of one geodesic band.

    The band points are projected into a local PCA plane, where a 2D circle is
    fitted. The circle center is then converted back into 3D.

    A residual-based median absolute deviation filter removes outliers and refits
    if enough points remain.
    """
    points = np.asarray(points, dtype=float)

    if len(points) < 3:
        raise RuntimeError("Too few points for section center.")

    # First fit.
    centroid, axis1, axis2, axis3 = fit_pca_basis(points)
    centered = points - centroid
    xy = np.column_stack([centered @ axis1, centered @ axis2])

    try:
        center_2d, radius_fit = fit_circle_2d(xy)
        center_3d = centroid + center_2d[0] * axis1 + center_2d[1] * axis2
    except Exception:
        # Fallback: use the band centroid if circle fitting fails.
        center_3d = centroid
        radius_fit = float(np.mean(np.linalg.norm(points - center_3d, axis=1)))

    distances = np.linalg.norm(points - center_3d, axis=1)

    # Robust outlier rejection using the median absolute deviation.
    median_distance = np.median(distances)
    mad = np.median(np.abs(distances - median_distance))

    if mad > 1e-12 and len(points) >= 20:
        keep = np.abs(distances - median_distance) <= 3.5 * 1.4826 * mad

        if np.sum(keep) >= 8:
            points_refit = points[keep]
            centroid, axis1, axis2, axis3 = fit_pca_basis(points_refit)
            centered = points_refit - centroid
            xy = np.column_stack([centered @ axis1, centered @ axis2])

            try:
                center_2d, radius_fit = fit_circle_2d(xy)
                center_3d = centroid + center_2d[0] * axis1 + center_2d[1] * axis2
            except Exception:
                center_3d = points_refit.mean(axis=0)
                radius_fit = float(np.mean(np.linalg.norm(points_refit - center_3d, axis=1)))

            distances = np.linalg.norm(points_refit - center_3d, axis=1)

    # Radius estimate.
    if RADIUS_AS_MEAN_DISTANCE:
        radius = float(np.mean(distances))
    else:
        radius = float(radius_fit)

    # Optional sanity checks.
    if MIN_REASONABLE_RADIUS is not None and radius < MIN_REASONABLE_RADIUS:
        raise RuntimeError("Radius below MIN_REASONABLE_RADIUS.")

    if MAX_REASONABLE_RADIUS is not None and radius > MAX_REASONABLE_RADIUS:
        raise RuntimeError("Radius above MAX_REASONABLE_RADIUS.")

    return center_3d, radius


# ============================================================
# CENTERLINE FROM GEODESIC BANDS
# ============================================================

def centerline_from_geodesic_bands(vertices, edges, loop0, loop1):
    """
    Compute the centerline using geodesic bands.

    This is the core algorithm. It assigns each vertex a normalized longitudinal
    coordinate s based on surface distances from the two rims. Vertices with
    similar s-values form a circumferential band.
    """
    graph = build_weighted_vertex_graph(vertices, edges)

    print("Computing geodesic distance from rim 1...")
    d0 = dijkstra_from_sources(graph, loop0)

    print("Computing geodesic distance from rim 2...")
    d1 = dijkstra_from_sources(graph, loop1)

    # Normalized longitudinal coordinate.
    denom = d0 + d1
    good = np.isfinite(denom) & (denom > 1e-12)

    s = np.full(len(vertices), np.nan, dtype=float)
    s[good] = d0[good] / denom[good]

    centerline = np.full((N_SECTIONS, 3), np.nan, dtype=float)
    radii = np.full(N_SECTIONS, np.nan, dtype=float)

    base_width = 1.0 / max(N_SECTIONS - 1, 1)
    targets = np.linspace(0.0, 1.0, N_SECTIONS)

    for i, target in enumerate(targets):
        selected = None

        # Try increasingly wider bands until enough vertices are found.
        for expansion in range(MAX_BIN_EXPANSIONS + 1):
            half_width = 0.5 * base_width * (BIN_EXPANSION_FACTOR ** expansion)

            if i == 0:
                mask = good & (s <= target + half_width)
            elif i == N_SECTIONS - 1:
                mask = good & (s >= target - half_width)
            else:
                mask = good & (np.abs(s - target) <= half_width)

            idx = np.where(mask)[0]

            if len(idx) >= MIN_POINTS_PER_SECTION:
                selected = idx
                break

        if selected is None:
            print(f"Section {i}: not enough vertices after expansion")
            continue

        points_in_band = vertices[selected]

        try:
            center, radius = robust_section_center_and_radius(points_in_band)
            centerline[i] = center
            radii[i] = radius
        except Exception as e:
            print(f"Section {i}: fit failed: {e}")

    # Force endpoints to be the centroids of the actual open rims.
    endpoint0 = vertices[loop0].mean(axis=0)
    endpoint1 = vertices[loop1].mean(axis=0)
    centerline[0] = endpoint0
    centerline[-1] = endpoint1

    # Fill missing centerline points by interpolation along section index.
    valid_centers = np.where(np.all(np.isfinite(centerline), axis=1))[0]

    if len(valid_centers) < 2:
        raise RuntimeError("Could not compute enough valid section centers.")

    for dim in range(3):
        y = centerline[:, dim]
        good_y = np.isfinite(y)
        centerline[:, dim] = np.interp(
            np.arange(N_SECTIONS),
            np.where(good_y)[0],
            y[good_y]
        )

    # Fill missing radii similarly.
    good_radii = np.isfinite(radii)

    if np.sum(good_radii) >= 2:
        radii = np.interp(
            np.arange(N_SECTIONS),
            np.where(good_radii)[0],
            radii[good_radii]
        )
    elif np.sum(good_radii) == 1:
        radii[:] = radii[good_radii][0]

    # Smooth the sequence of center points.
    centerline = smooth_centerline(
        centerline,
        passes=SMOOTHING_PASSES,
        strength=SMOOTHING_STRENGTH
    )

    # Re-anchor endpoints after smoothing.
    centerline[0] = endpoint0
    centerline[-1] = endpoint1

    return centerline, radii, s


def smooth_centerline(centerline, passes=8, strength=0.55):
    """Smooth interior centerline points while preserving endpoints."""
    out = np.asarray(centerline, dtype=float).copy()

    for _ in range(passes):
        new_out = out.copy()

        for i in range(1, len(out) - 1):
            neighbor_average = 0.5 * (out[i - 1] + out[i + 1])
            new_out[i] = (1.0 - strength) * out[i] + strength * neighbor_average

        new_out[0] = out[0]
        new_out[-1] = out[-1]
        out = new_out

    return out


# ============================================================
# ARC ANGLE
# ============================================================

def compute_arc_angle(centerline):
    """
    Estimate how much the tube arcs.

    The centerline is projected into its best-fit plane and a circle is fitted.
    Two angle estimates are reported:
        1. length_based_angle = centerline_length / fitted_arc_radius
        2. endpoint_angle = central angle between the endpoint radii
    """
    centroid, axis1, axis2, axis3 = fit_pca_basis(centerline)

    centered = centerline - centroid
    xy = np.column_stack([centered @ axis1, centered @ axis2])

    circle_center, circle_radius = fit_circle_2d(xy)

    angles = np.arctan2(
        xy[:, 1] - circle_center[1],
        xy[:, 0] - circle_center[0]
    )

    unwrapped = np.unwrap(angles)
    endpoint_angle = abs(unwrapped[-1] - unwrapped[0])

    length = polyline_length(centerline)
    length_based_angle = length / circle_radius

    return {
        "arc_fit_radius": float(circle_radius),
        "endpoint_angle_deg": float(math.degrees(endpoint_angle)),
        "length_based_angle_deg": float(math.degrees(length_based_angle)),
    }


# ============================================================
# BLENDER OUTPUT HELPERS
# ============================================================

def create_curve_polyline(points, name, bevel_depth=0.01):
    """Create a visible Blender curve object from a sequence of 3D points."""
    curve = bpy.data.curves.new(name, type="CURVE")
    curve.dimensions = "3D"
    curve.resolution_u = 2
    curve.bevel_depth = bevel_depth
    curve.bevel_resolution = 3

    spline = curve.splines.new("POLY")
    spline.points.add(len(points) - 1)

    for p, co in zip(spline.points, points):
        p.co = (float(co[0]), float(co[1]), float(co[2]), 1.0)

    obj = bpy.data.objects.new(name, curve)
    bpy.context.collection.objects.link(obj)
    return obj


def basis_from_tangent(tangent):
    """Build two perpendicular directions around a centerline tangent."""
    tangent = normalize(tangent)

    if abs(tangent[0]) < 0.8:
        temp = np.array([1.0, 0.0, 0.0])
    else:
        temp = np.array([0.0, 1.0, 0.0])

    u = normalize(np.cross(tangent, temp))
    v = normalize(np.cross(tangent, u))

    return u, v


def tangent_at(centerline, i):
    """Estimate the tangent direction at centerline point i."""
    if i == 0:
        return normalize(centerline[1] - centerline[0])

    if i == len(centerline) - 1:
        return normalize(centerline[-1] - centerline[-2])

    return normalize(centerline[i + 1] - centerline[i - 1])


def create_circle_curve(center, tangent, radius, name, bevel_depth=0.004, n_points=96, collection=None):
    """Draw a circular cross-section object for visualization."""
    u, v = basis_from_tangent(tangent)

    theta = np.linspace(0.0, 2.0 * math.pi, n_points + 1)

    points = (
        center[None, :]
        + radius * np.cos(theta)[:, None] * u[None, :]
        + radius * np.sin(theta)[:, None] * v[None, :]
    )

    curve = bpy.data.curves.new(name, type="CURVE")
    curve.dimensions = "3D"
    curve.resolution_u = 2
    curve.bevel_depth = bevel_depth
    curve.bevel_resolution = 2

    spline = curve.splines.new("POLY")
    spline.points.add(len(points) - 1)

    for p, co in zip(spline.points, points):
        p.co = (float(co[0]), float(co[1]), float(co[2]), 1.0)

    obj = bpy.data.objects.new(name, curve)

    if collection is None:
        bpy.context.collection.objects.link(obj)
    else:
        collection.objects.link(obj)

    return obj


def create_section_circle_objects(centerline, radii):
    """Create visual radius circles along the centerline."""
    collection = bpy.data.collections.new(SECTION_COLLECTION_NAME)
    bpy.context.scene.collection.children.link(collection)

    for i in range(len(centerline)):
        if not np.isfinite(radii[i]):
            continue

        create_circle_curve(
            center=centerline[i],
            tangent=tangent_at(centerline, i),
            radius=radii[i],
            name=f"geodesic_section_{i:03d}",
            bevel_depth=SECTION_CIRCLE_BEVEL_DEPTH,
            collection=collection
        )


def create_endpoint_marker(location, name, radius=0.05):
    """Create a small sphere marker at an endpoint centroid."""
    bpy.ops.mesh.primitive_uv_sphere_add(
        segments=16,
        ring_count=8,
        radius=radius,
        location=(float(location[0]), float(location[1]), float(location[2]))
    )

    obj = bpy.context.object
    obj.name = name
    return obj


# ============================================================
# CSV OUTPUT HELPERS
# ============================================================

def save_section_csv(centerline, radii):
    """
    Save one row per centerline section.

    Columns:
        section_index
        x, y, z
        radius
    """
    filepath = os.path.join(get_output_directory(), SECTION_CSV_FILENAME)

    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["section_index", "x", "y", "z", "radius"])

        for i, p in enumerate(centerline):
            writer.writerow([
                i,
                float(p[0]),
                float(p[1]),
                float(p[2]),
                float(radii[i]) if np.isfinite(radii[i]) else ""
            ])

    print(f"Saved section CSV to: {filepath}")
    return filepath


def save_summary_csv(summary):
    """
    Save the console-style summary metrics to a separate CSV file.

    This directly addresses the need to preserve the numerical console outputs
    in a file. The file has two columns:
        metric
        value
    """
    filepath = os.path.join(get_output_directory(), SUMMARY_CSV_FILENAME)

    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])

        for key, value in summary.items():
            writer.writerow([key, value])

    print(f"Saved summary metrics CSV to: {filepath}")
    return filepath


# ============================================================
# MAIN
# ============================================================

def main():
    """Run centerline extraction on the selected tube mesh."""
    obj = get_selected_mesh_object()

    print("\n========================================")
    print("Geodesic-band centerline for hollow open tube wall")
    print("Object:", obj.name)
    print("========================================")

    vertices, faces, edges = get_world_vertices_faces_edges(obj)

    print(f"Vertices: {len(vertices)}")
    print(f"Faces:    {len(faces)}")
    print(f"Edges:    {len(edges)}")

    # Detect the two open rim loops.
    loops = find_boundary_loops(faces)

    print(f"Boundary loops found: {len(loops)}")
    print("Boundary loop sizes:", [len(loop) for loop in loops])

    if len(loops) != 2:
        raise RuntimeError(
            "Expected exactly two boundary loops. "
            "Repair holes or seams until non-manifold selection shows only the two open rims."
        )

    loop0 = loops[0]
    loop1 = loops[1]

    endpoint0 = vertices[loop0].mean(axis=0)
    endpoint1 = vertices[loop1].mean(axis=0)

    # Add endpoint markers for visual confirmation.
    if CREATE_ENDPOINT_MARKERS:
        create_endpoint_marker(endpoint0, "tube_end_1_centroid", ENDPOINT_MARKER_RADIUS)
        create_endpoint_marker(endpoint1, "tube_end_2_centroid", ENDPOINT_MARKER_RADIUS)

    # Compute centerline and local radii.
    centerline, radii, s_coordinate = centerline_from_geodesic_bands(
        vertices,
        edges,
        loop0,
        loop1
    )

    length = polyline_length(centerline)

    valid_radii = radii[np.isfinite(radii)]

    if len(valid_radii) == 0:
        raise RuntimeError("No valid radii computed.")

    angle_info = compute_arc_angle(centerline)

    # Build summary dictionary first, then print and save it.
    summary = {
        "object_name": obj.name,
        "vertices": int(len(vertices)),
        "faces": int(len(faces)),
        "edges": int(len(edges)),
        "boundary_loops_found": int(len(loops)),
        "boundary_loop_1_vertices": int(len(loop0)),
        "boundary_loop_2_vertices": int(len(loop1)),
        "n_centerline_sections": int(len(centerline)),
        "valid_radii_count": int(len(valid_radii)),
        "centerline_length_blender_units": float(length),
        "average_radius_blender_units": float(np.mean(valid_radii)),
        "median_radius_blender_units": float(np.median(valid_radii)),
        "radius_std_dev_blender_units": float(np.std(valid_radii)),
        "minimum_radius_blender_units": float(np.min(valid_radii)),
        "maximum_radius_blender_units": float(np.max(valid_radii)),
        "fitted_arc_radius_blender_units": float(angle_info["arc_fit_radius"]),
        "angle_from_length_over_arc_radius_deg": float(angle_info["length_based_angle_deg"]),
        "endpoint_central_angle_deg": float(angle_info["endpoint_angle_deg"]),
        "endpoint0_x": float(endpoint0[0]),
        "endpoint0_y": float(endpoint0[1]),
        "endpoint0_z": float(endpoint0[2]),
        "endpoint1_x": float(endpoint1[0]),
        "endpoint1_y": float(endpoint1[1]),
        "endpoint1_z": float(endpoint1[2]),
    }

    # Create Blender visualization objects.
    if CREATE_CENTERLINE_OBJECT:
        create_curve_polyline(
            centerline,
            CENTERLINE_OBJECT_NAME,
            CENTERLINE_BEVEL_DEPTH
        )

    if CREATE_SECTION_CIRCLES:
        create_section_circle_objects(centerline, radii)

    # Save CSV outputs.
    section_csv_path = None
    summary_csv_path = None

    if SAVE_SECTION_CSV:
        section_csv_path = save_section_csv(centerline, radii)

    if SAVE_SUMMARY_CSV:
        summary_csv_path = save_summary_csv(summary)

    # Print measurements to the console.
    print("\n========================================")
    print("RESULTS")
    print("========================================")
    print(f"Centerline length:       {summary['centerline_length_blender_units']:.6f} Blender units")
    print(f"Average radius:          {summary['average_radius_blender_units']:.6f} Blender units")
    print(f"Median radius:           {summary['median_radius_blender_units']:.6f} Blender units")
    print(f"Radius std dev:          {summary['radius_std_dev_blender_units']:.6f} Blender units")
    print(f"Minimum radius:          {summary['minimum_radius_blender_units']:.6f} Blender units")
    print(f"Maximum radius:          {summary['maximum_radius_blender_units']:.6f} Blender units")
    print(f"Valid radii:             {summary['valid_radii_count']} / {summary['n_centerline_sections']}")

    print("\nArc angle estimate:")
    print(f"  Fitted arc radius:              {summary['fitted_arc_radius_blender_units']:.6f} Blender units")
    print(f"  Angle from length / arc radius: {summary['angle_from_length_over_arc_radius_deg']:.3f} degrees")
    print(f"  Endpoint central angle:         {summary['endpoint_central_angle_deg']:.3f} degrees")

    print("\nObjects created:")
    if CREATE_ENDPOINT_MARKERS:
        print("  tube_end_1_centroid")
        print("  tube_end_2_centroid")
    if CREATE_CENTERLINE_OBJECT:
        print(f"  {CENTERLINE_OBJECT_NAME}")
    if CREATE_SECTION_CIRCLES:
        print(f"  Collection: {SECTION_COLLECTION_NAME}")

    print("\nCSV files:")
    if section_csv_path is not None:
        print(f"  Section data:   {section_csv_path}")
    if summary_csv_path is not None:
        print(f"  Summary metrics: {summary_csv_path}")

    print("========================================\n")

    # Store results in a dictionary for interactive inspection in Blender.
    return {
        "centerline": centerline,
        "radii": radii,
        "length": length,
        "average_radius": float(np.mean(valid_radii)),
        "median_radius": float(np.median(valid_radii)),
        "radius_std": float(np.std(valid_radii)),
        "angle_info": angle_info,
        "endpoint0": endpoint0,
        "endpoint1": endpoint1,
        "s_coordinate": s_coordinate,
        "summary": summary,
        "section_csv_path": section_csv_path,
        "summary_csv_path": summary_csv_path,
    }


# Run the script.
results = main()
