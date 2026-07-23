#!/usr/bin/env python3
"""
genFigure1a.py
===================================================================
Build the annotated, branch-coloured phylogeny directly from the
IQ-TREE Newick file (exact branch lengths; no PDF reverse-engineering).

Pipeline (all deterministic; no generative component):
  1. Read the Newick tree.
  2. Re-root it between Protostomia and Deuterostomia -- i.e. at the
     midpoint of the branch subtending the vertebrate clade -- so the
     display matches the published rooting.
  3. Order each node's children to reproduce the agreed top-to-bottom
     tip order (finer control than generic ladderisation).
  4. Lay out a rectangular phylogram: x = cumulative branch length,
     tips evenly spaced in y, internal y = midpoint of children.
  5. Read per-tip trait values (pS) and estimate internal-node values
     by Brownian-motion ancestral-state reconstruction (a linear solve).
  6. Render to PDF: branches coloured on a continuous 'turbo' scale,
     italic Latin tip labels, right-side clade brackets, a left-margin
     colour bar, and a scale bar.

Dependencies:  pip install dendropy reportlab matplotlib numpy
Usage:         python3 make_tree_figure_from_newick.py

Inputs (edit paths below):
  IN_TREE : Newick tree file (.tree / .treefile / .contree)
  VALUES  : whitespace-delimited table  ->  Animal  pS  phylum
Output:
  OUT_PDF : finished figure
===================================================================
"""

import numpy as np
import dendropy
from reportlab.pdfgen import canvas
from reportlab.lib.colors import black
from matplotlib import cm
from matplotlib.colors import Normalize
from matplotlib import pyplot

# ----------------------------------------------------------------- config
IN_TREE = "iqtree_LG+C20+F+G_20260624044521-renamed.tree"
VALUES  = "pSvalues.txt"
OUT_PDF = "Figures/Figure1a.pdf"

CMAP_NAME = "turbo"
LINEWIDTH = 4.0
FS        = 16          # italic Latin tip labels
CLADEFS   = 18          # clade bracket labels
CLADEGAP  = 24          # gap between bracket bar and its label (pt)
SCALEFS   = 18
CBTICK    = 15
CBTITLE   = 18

# --- layout constants (these reproduce the original figure's proportions;
#     change freely for a different tree) ---
PT_PER_SUBST = 150.48   # horizontal scale (points per substitution/site)
ROOT_X       = 12.67    # x of the root node
TIP_Y_TOP    = 11.0     # y of the top-most tip
TIP_SPACING  = 32.611   # vertical gap between adjacent tips
STUB         = 6.67     # little branch drawn to the left of the root
MARGIN       = 100.0    # left margin that holds the colour bar

# The clade the tree is rooted on (root is placed at this branch's midpoint):
ROOT_CLADE = ["Little_Skate", "Tiger_Salamander", "White_Leghorn_Chicken",
              "Cow", "Dog", "Human", "House_mouse", "Brown_rat", "Rabbit"]

# Desired top-to-bottom tip order (children are ordered to realise this):
TARGET_ORDER = [
    "Longfin_Inshore_squid", "California_sea_hare", "Burgundy_snail",
    "Garden_snail", "C_elegans", "Red_swamp_Crayfish", "Fruit_fly",
    "Desert_Locust", "Cricket", "American_Cockroach", "Little_Skate",
    "Tiger_Salamander", "White_Leghorn_Chicken", "Cow", "Dog", "Human",
    "House_mouse", "Brown_rat", "Rabbit",
]

LATIN = {
    "Longfin_Inshore_squid": "Doryteuthis pealeii",
    "California_sea_hare":    "Aplysia californica",
    "Garden_snail":           "Cornu aspersum",
    "Burgundy_snail":         "Helix pomatia",
    "C_elegans":              "Caenorhabditis elegans",
    "Red_swamp_Crayfish":     "Procambarus clarkii",
    "Fruit_fly":              "Drosophila melanogaster",
    "American_Cockroach":     "Periplaneta americana",
    "Desert_Locust":          "Schistocerca gregaria",
    "Cricket":                "Gryllus bimaculatus",
    "Little_Skate":           "Leucoraja erinacea",
    "Tiger_Salamander":       "Ambystoma tigrinum",
    "White_Leghorn_Chicken":  "Gallus gallus domesticus",
    "Dog":                    "Canis lupus familiaris",
    "Cow":                    "Bos taurus",
    "Human":                  "Homo sapiens",
    "Rabbit":                 "Oryctolagus cuniculus",
    "Brown_rat":              "Rattus norvegicus domestica",
    "House_mouse":            "Mus musculus",
}

CLADES_INNER = [
    ("Mollusca", ["Longfin_Inshore_squid", "California_sea_hare",
                  "Burgundy_snail", "Garden_snail"]),
    ("Insecta",  ["Fruit_fly", "Desert_Locust", "Cricket",
                  "American_Cockroach"]),
    ("Mammalia", ["Cow", "Dog", "Human", "House_mouse", "Brown_rat", "Rabbit"]),
]
CLADES_OUTER = [
    ("Protostomia", ["Longfin_Inshore_squid", "California_sea_hare",
                     "Burgundy_snail", "Garden_snail", "C_elegans",
                     "Red_swamp_Crayfish", "Fruit_fly", "Desert_Locust",
                     "Cricket", "American_Cockroach"]),
    ("Deuterostomia", ["Little_Skate", "Tiger_Salamander",
                       "White_Leghorn_Chicken", "Cow", "Dog", "Human",
                       "House_mouse", "Brown_rat", "Rabbit"]),
]


# ============================================================ tree handling
class Node:
    def __init__(self):
        self.name = None          # leaf label, else None
        self.bl = 0.0             # branch length to parent (substitutions)
        self.kids = []
        self.parent = None
        self.x = self.y = 0.0
        self.ps = None


def load_and_root(path):
    """Read Newick, re-root at midpoint of the ROOT_CLADE stem branch."""
    t = dendropy.Tree.get(path=path, schema="newick",
                          preserve_underscores=True)
    t.is_rooted = True
    want = set(ROOT_CLADE)
    node = None
    for nd in t.preorder_node_iter():
        if not nd.is_leaf() and set(l.taxon.label for l in nd.leaf_iter()) == want:
            node = nd
            break
    if node is None:
        raise SystemExit("Could not locate the root clade in the tree.")
    L = node.edge.length
    t.reroot_at_edge(node.edge, length1=L / 2.0, length2=L / 2.0,
                     update_bipartitions=False)
    return t


def to_node(dnode):
    """Convert a dendropy subtree to our Node structure."""
    n = Node()
    n.bl = dnode.edge.length or 0.0
    if dnode.is_leaf():
        n.name = dnode.taxon.label
    else:
        for c in dnode.child_nodes():
            ch = to_node(c); ch.parent = n; n.kids.append(ch)
    return n


def all_nodes(n):
    yield n
    for k in n.kids:
        yield from all_nodes(k)


def leaves(n):
    return [n] if not n.kids else [l for k in n.kids for l in leaves(k)]


# ================================================= ordering + phylogram layout
def order_and_layout(root):
    rank = {name: i for i, name in enumerate(TARGET_ORDER)}

    def min_rank(n):
        return rank[n.name] if not n.kids else min(min_rank(k) for k in n.kids)

    def sort_kids(n):
        for k in n.kids:
            sort_kids(k)
        n.kids.sort(key=min_rank)
    sort_kids(root)

    tips = leaves(root)                       # now in target order
    for i, t in enumerate(tips):
        t.y = TIP_Y_TOP + i * TIP_SPACING

    def set_xy(n, x_parent_right):
        n.x = x_parent_right + n.bl * PT_PER_SUBST
        for k in n.kids:
            set_xy(k, n.x)
        if n.kids:
            ys = [k.y for k in n.kids]
            n.y = (min(ys) + max(ys)) / 2.0
    root.x = ROOT_X
    for k in root.kids:
        set_xy(k, root.x)
    if root.kids:
        ys = [k.y for k in root.kids]
        root.y = (min(ys) + max(ys)) / 2.0
    return tips


# ============================================ Brownian-motion ancestral states
def brownian_asr(root):
    internal = [n for n in all_nodes(root) if n.kids]
    idx = {id(n): i for i, n in enumerate(internal)}
    m = len(internal)
    A, b = np.zeros((m, m)), np.zeros(m)

    def w(n):
        return 1.0 / max(n.bl, 1e-6)
    for n in internal:
        i = idx[id(n)]; tot = 0.0
        nbrs = ([(n.parent, w(n))] if n.parent is not None else []) + \
               [(c, w(c)) for c in n.kids]
        for nb, ww in nbrs:
            tot += ww
            if nb.kids:
                A[i, idx[id(nb)]] -= ww
            else:
                b[i] += ww * nb.ps
        A[i, i] += tot
    z = np.linalg.solve(A, b)
    for n in internal:
        n.ps = float(z[idx[id(n)]])


# ============================================================= rendering
def render(root, tips, page_h, out_pdf):
    cmap = pyplot.get_cmap(CMAP_NAME)
    vmin, vmax = min(t.ps for t in tips), max(t.ps for t in tips)
    norm = Normalize(vmin, vmax)

    def color(v):
        r, g, b, _ = cmap(norm(v)); return (r, g, b)

    def X(x):
        return x + MARGIN

    off = FS * 0.33                                  # vertical label centring

    meas = canvas.Canvas("._m.pdf")
    maxend = max(t.x + 7 + meas.stringWidth(LATIN[t.name], "Helvetica-Oblique",
                 FS) for t in tips)
    inner_bx = maxend + 26;  inner_lx = inner_bx + CLADEGAP
    outer_bx = inner_bx + 44; outer_lx = outer_bx + CLADEGAP
    page_w = X(outer_lx) + 26

    c = canvas.Canvas(out_pdf, pagesize=(page_w, page_h)); c.setLineCap(1)

    def Y(v):
        return page_h - v

    def grad(xL, xR, y, vL, vR, n=44):
        for i in range(n):
            ta, tb = i / n, (i + 1) / n
            xa, xb = X(xL + ta * (xR - xL)), X(xL + tb * (xR - xL))
            r, g, b = color(vL + (ta + tb) / 2 * (vR - vL))
            c.setStrokeColorRGB(r, g, b); c.line(xa, Y(y), xb, Y(y))

    nodes = list(all_nodes(root))
    c.setLineWidth(LINEWIDTH)
    r, g, b = color(root.ps); c.setStrokeColorRGB(r, g, b)
    c.line(X(root.x - STUB), Y(root.y), X(root.x), Y(root.y))     # root stub
    for n in nodes:
        if n.parent is not None:
            grad(n.parent.x, n.x, n.y, n.parent.ps, n.ps)
        if n.kids:
            ys = [k.y for k in n.kids]
            r, g, b = color(n.ps); c.setStrokeColorRGB(r, g, b)
            c.line(X(n.x), Y(min(ys)), X(n.x), Y(max(ys)))

    # scale bar (0.5 substitutions/site) under the tree
    sb_len = 0.5 * PT_PER_SUBST
    sb_x = X(ROOT_X + 4.2 * PT_PER_SUBST)
    sb_y = TIP_Y_TOP + (len(tips) - 0.6) * TIP_SPACING
    c.setLineWidth(1.0); c.setStrokeColor(black)
    c.line(sb_x, Y(sb_y), sb_x + sb_len, Y(sb_y))
    c.setFont("Helvetica", SCALEFS); c.setFillColor(black)
    c.drawCentredString(sb_x + sb_len / 2, Y(sb_y + SCALEFS + 2), "0.5")

    # italic Latin tip labels
    c.setFont("Helvetica-Oblique", FS); c.setFillColor(black)
    for t in tips:
        c.drawString(X(t.x + 7), Y(t.y + off), LATIN[t.name])

    # clade brackets
    def bracket(bx, lx, name, members):
        ys = [tt.y for tt in tips if tt.name in members]
        y0, y1 = min(ys), max(ys)
        c.setStrokeColor(black); c.setLineWidth(1.4)
        c.line(X(bx), Y(y0), X(bx), Y(y1))
        c.line(X(bx), Y(y0), X(bx) - 6, Y(y0))
        c.line(X(bx), Y(y1), X(bx) - 6, Y(y1))
        c.saveState(); c.translate(X(lx), Y((y0 + y1) / 2)); c.rotate(90)
        c.setFont("Helvetica-Bold", CLADEFS); c.setFillColor(black)
        c.drawCentredString(0, 0, name); c.restoreState()
    for nm, mem in CLADES_INNER:
        bracket(inner_bx, inner_lx, nm, mem)
    for nm, mem in CLADES_OUTER:
        bracket(outer_bx, outer_lx, nm, mem)

    # colour bar in the LEFT margin, vertically centred, ticks on the right
    cbx, cby, h, w, M = 24, 213, 215, 16, 150
    for i in range(M):
        v = vmin + i / (M - 1) * (vmax - vmin)
        r, g, b = color(v); c.setFillColorRGB(r, g, b)
        c.rect(cbx, Y(cby + h - i * (h / M)), w, h / M + 0.8, stroke=0, fill=1)
    c.setStrokeColor(black); c.setFillColor(black); c.setLineWidth(0.8)
    c.rect(cbx, Y(cby + h), w, h, stroke=1, fill=0)
    c.setFont("Helvetica", CBTICK)
    for v, lab in [(40, "40"), (80, "80"), (120, "120"),
                   (160, "160"), (200, "200"), (vmax, "254")]:
        fr = (v - vmin) / (vmax - vmin); yy = cby + h - fr * h
        c.line(cbx + w, Y(yy), cbx + w + 3, Y(yy))
        c.drawString(cbx + w + 6, Y(yy + 4.5), lab)
    c.setFont("Helvetica-Bold", CBTITLE)
    c.drawString(cbx - 2, Y(cby - 11), "pS")

    c.showPage(); c.save()
    print(f"wrote {out_pdf}  ({page_w:.0f} x {page_h:.0f} pt)")


# ==================================================================== main
def main():
    dtree = load_and_root(IN_TREE)
    root = to_node(dtree.seed_node)
    tips = order_and_layout(root)

    vals = {}
    for ln in open(VALUES):
        p = ln.split()
        if not p or p[0] == "Animal":
            continue
        vals[p[0]] = float(p[1])
    for t in tips:
        t.ps = vals[t.name]

    brownian_asr(root)

    page_h = TIP_Y_TOP + (len(tips) - 1) * TIP_SPACING + 44   # bottom padding
    render(root, tips, page_h, out_pdf=OUT_PDF)


if __name__ == "__main__":
    main()
