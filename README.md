<p align="center">
  <img src="utilities/logo/convert3d.jpg" alt="Convert3D" width="600">
</p>

# Convert3D

Convert3D (`c3d`) is a command-line image processing tool that complements
[ITK-SNAP](http://www.itksnap.org). It was originally developed to convert 3D medical
images between file formats. It has since become a Swiss army knife for medical image
processing, with commands for resampling, smoothing, thresholding, computing features
and image statistics, and more. It also has commands for multi-label images, such as
ITK-SNAP segmentations, and for multi-component images, such as RGB images.

c3d works like a [Reverse Polish notation](https://en.wikipedia.org/wiki/Reverse_Polish_notation)
calculator. Images are placed on a stack, and each command acts on the images at the top
of the stack. This lets you chain several operations into a single command without
saving intermediate files. For example, this command smooths an image, resamples it to
half its resolution and saves the result:

```
c3d input.nii.gz -smooth 2mm -resample 50% -o output.nii.gz
```

The companion tools `c2d` and `c4d` provide the same commands for 2D and 4D images.

## Getting Convert3D

- **ITK-SNAP:** c3d comes with [ITK-SNAP](http://www.itksnap.org). Choose
  *Help → Install Command-Line Tools* from the ITK-SNAP menu.
- **Binaries:** download them from [SourceForge](https://sourceforge.net/projects/c3d/).
- **Source:** build with CMake against [ITK](https://itk.org).
- **Python:** `pip install picsl_c3d` installs the Python wrapper
  ([c3d_python](https://github.com/pyushkevich/c3d_python)). It runs c3d commands from
  Python and exposes the images on the c3d stack as [SimpleITK](https://simpleitk.org)
  images.

## Documentation

- [Reference manual](doc/c3d.md), covering every command and option.
- [Convert3D page](https://www.itksnap.org/pmwiki/pmwiki.php?n=Convert3D.Convert3D) on
  itksnap.org.
- [README.txt](README.txt) explains how to add a new command.

The logo is made from real c3d outputs; see [utilities/logo](utilities/logo).

Convert3D is developed at the [Penn Image Computing and Science Laboratory
(PICSL)](https://picsl.upenn.edu) at the University of Pennsylvania.
