"""Direct labels that do not overprint each other or the data.

Points in these figures crowd (a log-distributed community along a diagonal,
a dozen organisms near Gini = 1), so a fixed offset, or alternating sides, still
lets neighbours collide. place_labels() tries a ring of positions around each
point, nearest first, then two rings further out with a thin leader line, and
keeps the first position whose text box clears every label already placed,
every marker, any extra obstacle (error bars, connecting lines) and any artist
named in `avoid` (a legend, a note), and stays inside the axes. If nothing
clears, it takes the position with the least overlap.

aggregate_results.py carries its own copy of this function, because the
AGGREGATE process stages that script alone; keep the two in step.
"""
from __future__ import annotations


def place_labels(ax, points, texts, fontsize=6.5, color="0.25", style="normal",
                 obstacles=(), avoid=(), marker_px=5):
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    to_px = ax.transData.transform
    axbox = ax.get_window_extent(renderer)
    taken = list(obstacles)
    for artist in avoid:
        bb = artist.get_window_extent(renderer)
        taken.append((bb.x0, bb.y0, bb.x1, bb.y1))
    for (x, y) in points:
        px, py = to_px((x, y))
        taken.append((px - marker_px, py - marker_px, px + marker_px, py + marker_px))

    def overlap(a, b):
        w = min(a[2], b[2]) - max(a[0], b[0])
        h = min(a[3], b[3]) - max(a[1], b[1])
        return w * h if w > 0 and h > 0 else 0.0

    base = [(6, 4, "left", "bottom"), (6, -4, "left", "top"),
            (-6, 4, "right", "bottom"), (-6, -4, "right", "top"),
            (8, 0, "left", "center"), (-8, 0, "right", "center"),
            (0, 8, "center", "bottom"), (0, -8, "center", "top"),
            (4, 10, "left", "bottom"), (4, -10, "left", "top"),
            (-4, 10, "right", "bottom"), (-4, -10, "right", "top")]
    cands = [(dx * k, dy * k, ha, va, k > 1)
             for k in (1, 2.2, 3.5) for dx, dy, ha, va in base]
    pad = 1.5 * fig.dpi / 72.0
    placed = []
    for (x, y), text in zip(points, texts):
        best = None
        for dx, dy, ha, va, leader in cands:
            # Measured without the leader line: an annotation's box includes
            # its arrow, which reaches back to the point and would count the
            # label's own marker as an overlap.
            ann = ax.annotate(text, (x, y), textcoords="offset points",
                              xytext=(dx, dy), ha=ha, va=va, fontsize=fontsize,
                              color=color, style=style)
            bb = ann.get_window_extent(renderer)
            ann.remove()
            box = (bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad)
            outside = (box[0] < axbox.x0 or box[2] > axbox.x1
                       or box[1] < axbox.y0 or box[3] > axbox.y1)
            cost = sum(overlap(box, t) for t in taken) + (1e9 if outside else 0)
            if best is None or cost < best[0]:
                best = (cost, (dx, dy, ha, va, leader), box)
            if cost == 0:
                break
        dx, dy, ha, va, leader = best[1]
        placed.append(ax.annotate(
            text, (x, y), textcoords="offset points", xytext=(dx, dy), ha=ha, va=va,
            fontsize=fontsize, color=color, style=style,
            arrowprops=(dict(arrowstyle="-", color="0.6", lw=0.5, shrinkA=1, shrinkB=4)
                        if leader else None)))
        taken.append(best[2])
    return placed


def segment_boxes(ax, p0, p1, half_px=2, step_px=6):
    """Small pixel boxes along a segment between two data points.

    A chain of small boxes, not one bounding box: a sloped segment's bounding
    box would block a whole rectangle the line never crosses.
    """
    import math
    (x0, y0), (x1, y1) = ax.transData.transform(p0), ax.transData.transform(p1)
    n = max(1, int(math.hypot(x1 - x0, y1 - y0) / step_px))
    out = []
    for i in range(n + 1):
        x, y = x0 + (x1 - x0) * i / n, y0 + (y1 - y0) * i / n
        out.append((x - half_px, y - half_px, x + half_px, y + half_px))
    return out
