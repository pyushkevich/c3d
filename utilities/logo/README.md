# Convert3D logo

`convert3d.jpg` is the Convert3D logo, and `convert3d-card.jpg` is a smaller version
without captions for places where captions would be too small to read, such as website
cards. Both are made from real c3d outputs. A crop of a T1 MRI is drawn as a 3D box,
next to boxes showing the outputs of `-resample 20%`, `-smooth 2mm`, `-grad`, `-otsu 3`
and `-sdt`.

To regenerate them, you need `c3d` on the path, Python 3 with `numpy`, `nibabel` and
`pillow`, the MRI
[`MRIcrop-orig.gipl.gz`](https://github.com/pyushkevich/itksnap/raw/master/Testing/TestData/MRIcrop-orig.gipl.gz)
from the ITK-SNAP test data, and the font
[Barlow SemiBold](https://github.com/google/fonts/raw/main/ofl/barlow/Barlow-SemiBold.ttf)
(SIL Open Font License). Then run:

```
./make_logo.sh MRIcrop-orig.gipl.gz work Barlow-SemiBold.ttf
```

This writes the c3d outputs and both images to `work/`.

- `make_logo.sh` runs the c3d commands and then `compose_logo.py`.
- `compose_logo.py` draws the boxes.
