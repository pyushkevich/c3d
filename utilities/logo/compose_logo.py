"""Compose a Convert3D logo from real c3d outputs (see make_logo.sh), drawing each
output as a cut-away box showing its axial, coronal and sagittal faces."""
import sys
import numpy as np
import nibabel as nib
from PIL import Image, ImageDraw, ImageFont, ImageFilter

src, font, out = sys.argv[1], sys.argv[2], sys.argv[3]
THUMB = len(sys.argv) > 4 and sys.argv[4] == "thumb"  # small version without captions
DARK, LIGHT = (0x16, 0x56, 0x7a), (0x00, 0xa8, 0xd8)
GOLD = (0xf4, 0xc4, 0x5a)
TEXT = (0x33, 0x33, 0x33)

# box cut from the 78 x 110 x 64 volume; top face is the axial slice at z = Z1
X0, X1, Y0, Y1, Z0, Z1 = 15, 63, 31, 49, 2, 64


def load(name):
    return nib.load(f"{src}/{name}.nii.gz").get_fdata().squeeze()


def grey(lo=0.5, hi=99.6):
    def f(v):
        a, b = np.percentile(v, [lo, hi])
        g = np.clip((v - a) / (b - a), 0, 1) * 255
        return np.stack([g] * 3, -1)
    return f


def labels(v):
    cols = np.array([(10, 40, 60), DARK, LIGHT, GOLD], float)
    return cols[v.astype(int)]


def distance(v):
    t = np.clip(v / 18, 0, 1)[..., None]
    rgb = np.array(DARK) + (np.array((235, 248, 252)) - np.array(DARK)) * t
    rgb = np.where((v < 0)[..., None], np.array(GOLD, float), rgb)
    bands = (np.abs(np.mod(v, 3) - 1.5) < 0.2) & (v > 0)
    rgb[bands] *= 0.75
    return rgb


def faces(vol, colorize):
    rgb = colorize(vol[X0:X1, Y0:Y1, Z0:Z1])
    front = rgb[:, 0, ::-1].transpose(1, 0, 2)    # coronal face j = 0: rows top -> bottom
    top = rgb[:, ::-1, -1].transpose(1, 0, 2)     # rows back -> front, cols i
    side = rgb[-1, :, ::-1].transpose(1, 0, 2)    # face i = nx: rows top -> bottom, cols j
    return [f.astype(np.uint8) for f in (front, top, side)]


def warp(face, origin, u, v, size, nearest):
    """Map face columns along vector u and rows along v, starting at origin."""
    h, w = face.shape[:2]
    im = Image.fromarray(face)
    M = np.array([[u[0], v[0]], [u[1], v[1]]], float)
    Mi = np.linalg.inv(M)
    # output (X, Y) -> (s, t) in [0,1]^2 -> input pixel (s*w, t*h)
    a, b = Mi[0] * w
    d, e = Mi[1] * h
    c = -(a * origin[0] + b * origin[1])
    f = -(d * origin[0] + e * origin[1])
    res = Image.NEAREST if nearest else Image.BICUBIC
    layer = im.transform(size, Image.AFFINE, (a, b, c, d, e, f), resample=res)
    mask = Image.new("L", size, 0)
    p = np.array(origin, float)
    ImageDraw.Draw(mask).polygon([tuple(p), tuple(p + u), tuple(p + u + v), tuple(p + v)], fill=255)
    return layer, mask


DEPTH = np.array([0.42, -0.3])  # screen direction of the j (depth) axis


def box(vol, colorize, o, scale):
    """Upright slab: coronal face toward the viewer with its top left corner at o."""
    front, top, side = faces(vol, colorize)
    nx, ny, nz = X1 - X0, Y1 - Y0, Z1 - Z0
    a = np.array([1.0, 0]) * scale * nx
    b = DEPTH * scale * ny
    c = np.array([0, 1.0]) * scale * nz
    o = np.array(o, float)
    return [
        (top, o + b, a, -b, 0.82),
        (side, o + a, b, c, 0.68),
        (front, o, a, c, 1.0),
    ]


W, H = 1000, 500
canvas = Image.new("RGB", (W, H), "white")
draw = ImageDraw.Draw(canvas)
def mono_font(size):
    for path in ("/System/Library/Fonts/Menlo.ttc", "DejaVuSansMono.ttf", "consola.ttf"):
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            pass
    return ImageFont.load_default(size)


mono = mono_font(13)


def render(vol, colorize, o, scale, nearest=False, lean=0.0):
    """Draw a slab, tilted clockwise by `lean` degrees about the bottom right
    corner of its front face, like a book leaning on its neighbour."""
    quads = box(vol, colorize, o, scale)
    layer = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    ld = ImageDraw.Draw(layer)
    for face, org, u, v, shade in quads:
        im, mask = warp((face * shade).astype(np.uint8), org, u, v, canvas.size, nearest)
        layer.paste(im, (0, 0), mask)
    for _, org, u, v, _ in quads:
        q = [org, org + u, org + u + v, org + v, org]
        ld.line([tuple(p) for p in q], fill=(34, 34, 34, 255), width=1, joint="curve")
    _, org, a, c, _ = quads[2]
    pivot = tuple(org + a + c)
    layer = layer.rotate(-lean, resample=Image.BICUBIC, center=pivot)
    alpha = layer.getchannel("A")
    shadow = Image.new("L", canvas.size, 0)
    shadow.paste(alpha, (5, 8))
    shadow = shadow.filter(ImageFilter.GaussianBlur(9)).point(lambda x: x * 0.3)
    canvas.paste(Image.new("RGB", canvas.size, (20, 40, 60)), (0, 0), shadow)
    canvas.paste(layer, (0, 0), layer)


captions = []


def caption(text, cx, y):
    # drawn at the end; the thumbnail leaves them out (they would be too small to read)
    captions.append((text, cx, y))


def arrow(x, y, length):
    """Draw an arrow at 4x resolution and downsample it, since PIL does not antialias."""
    k, pad = 4, 12
    x0, y0 = int(x) - pad, int(y) - pad
    w, h = int(length) + 2 * pad, 2 * pad
    layer = Image.new("RGBA", (w * k, h * k), LIGHT + (0,))
    d = ImageDraw.Draw(layer)
    X, Y = (x - x0) * k, (y - y0) * k
    L = length * k
    d.rectangle((X + 4 * k, Y - 2.5 * k, X + L - 14 * k, Y + 2.5 * k), fill=LIGHT + (255,))
    d.polygon([(X + L - 4 * k, Y), (X + L - 18 * k, Y - 10 * k), (X + L - 18 * k, Y + 10 * k)], fill=LIGHT + (255,))
    layer = layer.resize((w, h), Image.LANCZOS)
    canvas.paste(layer, (x0, y0), layer)


nx, ny, nz = X1 - X0, Y1 - Y0, Z1 - Z0


def extent(scale):
    # bounding box relative to the top left corner of the front face
    return 0, (nx + DEPTH[0] * ny) * scale, DEPTH[1] * ny * scale, nz * scale


panels = [
    ("-resample 20%", "resample", grey(), True),
    ("-smooth 2mm", "smooth", grey(), False),
    ("-grad", "grad", grey(1, 99), False),
    ("-otsu 3", "otsu", labels, True),
    ("-sdt", "sdt", distance, False),
]

# six slabs of the same size standing on a shelf: the input, an arrow, then the
# outputs, the first two leaning on their neighbours like books
s, gap, arrow_w = 2.2, 30, 50  # gap includes the side face
leans = [0, 12, 6, 0, 0, 0]
l, r, t, b = extent(s)
hf, wf = nz * s, nx * s              # front face size
lean_dx = [hf * np.sin(np.radians(q)) for q in leans]
widths = [wf + max(0, r - wf), *[wf + max(gap, d + 10) for d in lean_dx[1:-1]], r]
total = widths[0] + arrow_w + sum(widths[1:])
x = (W - total) / 2
shelf = 270                          # y of the front faces' bottom edge
boxes = [("input.nii.gz", "original", grey(), False)] + panels
x_first = x
for i, (text, name, col, nearest) in enumerate(boxes):
    render(load(name), col, (x, shelf - hf), s, nearest, leans[i])
    caption(text, x + wf / 2, shelf + 14)
    if i == 0:
        x += widths[0]
        arrow(x, shelf - hf / 2, arrow_w)
        x += arrow_w
    else:
        x += widths[i] if i < len(widths) else 0
right = x
draw.rounded_rectangle((x_first - 14, shelf + 1, right + 14, shelf + 6), 3, fill=(214, 224, 231))

f = ImageFont.truetype(font, 132)
w1, w2 = draw.textlength("Convert", font=f), draw.textlength("3D", font=f)
tx = (W - w1 - w2) / 2
draw.text((tx, 300), "Convert", font=f, fill=DARK)
draw.text((tx + w1, 300), "3D", font=f, fill=LIGHT)
if THUMB:
    # drop the captions (too small to read) and the space they take, stack the boxes
    # tightly over the wordmark, then pad to 2:1 and size for website cards (330 x 165
    # px on the page, 2x for high DPI screens)
    def bbox(top, bottom):
        a = np.asarray(canvas.convert("L"))[top:bottom] < 248
        rows, cols = np.nonzero(a.any(1))[0], np.nonzero(a.any(0))[0]
        return cols[0], top + rows[0], cols[-1] + 1, top + rows[-1] + 1
    shelves, mark = canvas.crop(bbox(0, shelf + 7)), canvas.crop(bbox(shelf + 40, H))
    gap, margin = 22, 24
    cw = max(shelves.width, mark.width) + 2 * margin
    ch = shelves.height + gap + mark.height + 2 * margin
    cw, ch = max(cw, 2 * ch), max(ch, cw // 2)
    thumb = Image.new("RGB", (cw, ch), "white")
    y0 = (ch - shelves.height - gap - mark.height) // 2
    thumb.paste(shelves, ((cw - shelves.width) // 2, y0))
    thumb.paste(mark, ((cw - mark.width) // 2, y0 + shelves.height + gap))
    thumb.resize((660, 330), Image.LANCZOS).save(out, quality=92)
    sys.exit()
for text, cx, y in captions:
    w = draw.textlength(text, font=mono)
    draw.text((cx - w / 2, y), text, font=mono, fill=TEXT)
canvas.save(out, quality=92)
