# [Convert3D][1] Documentation

[TOC]

### What's New?

*   **-cos**,**-sin**,**-atan2** commands 
*   **-align-landmarks** command 
*   **-color-map** command 
*   **-min**, **-max** commands 
*   **-accum**, **-endaccum** loop structure 
*   **-tile** command, great for stacking TIFFs into a 3D volume 
*   **-test-xxx** set of commands 
*   **-holefill** command 
*   New **c4d** command 
*   **-slice** command extended to support range of slices (e.g., *5:10*, *0:-1*, *0:2:-1*) 

### About convert3d 

**Convert3d** is a command-line tool for converting 3D images between common file formats. The tool also includes a growing list of commands for image manipulation, such as thresholding and resampling. The tool can also be used to obtain information about image files. 

The simplest way to use **convert3d** is to convert images between formats. Here are some examples. Note that you must use the `-o` command to specify the output of the conversion. 

    c3d input.img -o output.vtk
    c3d hello.cub -o hello.img
    c3d float.img -type short -o short.img

**Convert3d** works in a way like a [Revese Polish notation][2] calculator. This means that actions must be preceeded by their inputs: instead of writing 'a + b', we must write 'a b +'. In our case, the inputs are image files, and the actions are specified on the command line. For example, to add a pair of images, we would issue the following command: 

    c3d input1.img input2.img -add -o output.img

Here, `input1.img` and `input2.img` are a pair of image files (in Analyze format), while `-add` and `-o` are commands. The command `-add` takes two input images, adds them, and generates a single output image. The command `-o` (short for output) takes one input image and writes it to a filename specified after the command (in this case, `output.img`). 

The behavior of **convert3d** commands can be affected by options. For example, the `-resample` command is affected by the `-interpolation` option. In order for the option to affect a command, *the option must precede the command on the command line*. For example: 

    c3d input.img -interpolation Cubic -resample 50x30x40vox -o output.img

All image processing is done in double-precision reals. However you can use the **-type** option to save as any data type. By default, image intensities are rounded to the nearest integer (not rounded down) when saving as a integral data type. 

### Command autocompletion for Bash users

Users of the bash shell can take advantage of its command autocompletion features. When this is enabled, you can type on the command line 

    c3d input.img -re[TAB]

and as you hit the TAB key, the list of possible completions will be presented to you. This way you will have to go back to this reference less often. To enable this feature, download the script [bashcomp.sh][3] and add the following line to your **.profile** script; 

    source SOMEDIR/bashcomp.sh

where SOMEDIR should be replaced by the directory where you saved the script. 

### Coordinate Frames and NIFTI

**Convert3d** is NIFTI-compatible. The transformation between the image coordinates and physical coordinates is read from the NIFTI header and maintained through the pipeline. The physical coordinates are in the **RAS** frame: *x* increases from left to right, *y* from posterior to anterior, *z* from inferior to superior. You can display the coordinate transformation using the **-info-full** command. 

    $ c3d.exe input.nii -info-full
    Image #1:
      Image Dimensions   : [176, 255, 216]
      Bounding Box       : {[87.2353 103.015 -108.064], [263.235 358.015 107.936]}
      Voxel Spacing      : [1, 1, 1]
      Intensity Range    : [5, 1775]
      Mean Intensity     : 103.175
      Voxel->RAS x-form  : 
             1.00000     -0.00000     -0.00000    -87.23527 
            -0.00000      1.00000     -0.00000   -103.01535 
             0.00000      0.00000      1.00000   -108.06398 
             0.00000      0.00000      0.00000      1.00000 
      Image Metadata: 
        ITK_InputFilterName = NiftiImageIO
        ITK_OnDiskStorageTypeName = short

### Usage Examples

To convert an image from Metaimage format to short-valued Analyze format use 

    c3d input.mha -type short -o output.img.gz

To apply some arithmetic to the image, do the following 

    c3d input.mha -scale 4096 -shift 2048 -type short -o output.img.gz

To resample and Gaussian blur an image the following command is used 

    c3d input.hdr -interpolation Linear -resample 256x256x192 
        -smooth 3vox -o output.img.gz

To compute **a+2b** with a pair of images, do 

    c3d a.img b.img -scale 2 -add -o c.img

### Multilabel Images

Multilabel (ML) images are images where each voxel is assigned a label, chosen from a relatively small set of labels. Labels are typically integers, although they can be floating point numbers. For example, labels can be used to identify anatomical structures in an image: label 1 may correspond to the white matter, 2 to the gray matter, 3 to CSF, and 0 to the background. 

When processing ML images, one must take certain care to make sure that numerical operations do not blend labels into meaningless values. For instance, smoothing a ML image will create intermediate values at places where two labels are near each other. A similar issue comes up when resampling or reslicing these images. However, **convert3d** provides commands that split the ML image into a set of binary images (one for each label), allowing you to perform operations on each label separately. Subsequently, the processed label-wise images can be combined into a new ML image. This is accomplished using commands **-split** and **-merge**. In addition, the **-foreach** ... **-endfor** command structure allows operations to be applied sequentially to all images on the stack. 

For example, to smooth the components of a ML image separately, we can use the following command 

    c3d ml.img -split -foreach -smooth 3mm -endfor -merge -o mlsmooth.img

The **-merge** command implements voting between individual label images. A label assigned to a voxel during merge is that of the label image with highest intensity at that voxel. The values of the labels assigned during the **-merge** command are remembered from the time the **-split** command is run. So if the image *ml.img* has labels *1, 4, 10* then the output image will also have the same labels. 

When applying spatial transformations to ML images, one may choose to use nearest neighbor interpolation: 

    c3d ref.img ml.img -interpolation NearestNeighbor -reslice-matrix affine.mat -o mltran.img

However, this causes some aliasing of the results. It is sometimes useful to assume that labels are *fuzzy*, and to apply linear or cubic interpolation to the fuzzy labels. Here is how we can apply an affine transformation to a ML image this way: 

    c3d ref.img -popas ref ml.img -split -foreach -smooth 3mm 
      -insert ref 1 -reslice-matrix affine.mat -endfor -merge-o mltran.img

This command is somewhat complex, mainly because the **-reslice-matrix** command requires a reference image as the first operand, and we have to use named images to insert it on the stack in the right place. 

### Commands: Image Input/Output and Information

Passing an image on the command line (without any switches) results in that image being read and placed at the end of the stack. For example:

    c3d myimage.nii -info another.nii -info

will result in information printed for both images. At the end, `myimage.nii` will be in the first position on the stack and `another.nii` will be at the end of the stack.

#### -dicom-series-list: List image series in a DICOM directory

Syntax: `-dicom-series-list <directory>

Prints out a table of DICOM series ids and corresponding image information to standard output.

#### -dicom-series-read: Read a DICOM image series

Syntax: `-dicom-series-read <directory> <series_id>`

Imports a specific DICOM image series from a directory containing DICOM files. The **directory** parameter may also point to one of the DICOM files in the directory.
The **seried_id** is a string identifier for the series that can be obtained by calling **-dicom-series-list**

#### -info: Display brief image information        

Syntax: `-info`

Prints brief information about the last image on the stack. Does not affect the stack.

    c3d image.hdr -info

Use with the **-foreach** command to get information on multiple images

    c3d images*.nii -foreach -info -endfor

#### -info-full: Display verbose image information        

Syntax: `-info-full`

Prints extended information about the last image on the stack, such as the metadata dictionary. For example, 

    c3d image.hdr -info-full

#### -mcs, -multicomponent-split: Enable splitting of multi-component images on read

Syntax: `-mcs`

Enable reading of multi-component images. By default, when a multi-component image is encountered, the components are combined into a single image. Setting the **-mcs** flag changes this behavior, and each of the components is loaded sequentially. See the section below on multi-component image support. 

    $ c3d -mcs rgb.mha -foreach -probe 110x110x80mm -endfor
    Interpolated image value at 110 110 80 is 1
    Interpolated image value at 110 110 80 is 66
    Interpolated image value at 110 110 80 is 29

    $ c3d rgb.mha -foreach -probe 110x110x80mm -endfor
    Interpolated image value at 110 110 80 is 49.5198

#### -nomcs, -no-multicomponent-split: Disable splitting of multi-component images on read

Syntax: `-nomcs`

Used to reverse the effect of previous **-mcs** command.

#### -o: Output (write) last image on the stack to image file    

Syntax: `-o filename`

Write image, overriding an existing image. Without the **-o** option, **convert3d** will write an image only if it does not exist. The **-o** options protects input images from being accidentally deleted. Here we copy an image, changing format:

    c3d image1.mha -o image2.nii

The **-o** option can also be used to save an intermediate image in the stack: 

    c3d image1.img -threshold 1 10 1 0 -o thresh.img -resample 50% -o final.img


#### -omc, -output-multicomponent: Output multiple images to single file

Syntax: `-omc [number] filename`

Write multiple images on the **Convert3d** stack as a single multi-component image file. If the optional number *n* is specified, only the last *n* images on the stack will be used. Not all file formats support multi-component output. NIFTI is the safest bet.

    c3d red.nii green.nii blue.nii -omc rgb.mha

For 2D images, this command can be used to generate color PNG files:

    c3d image.nii -slice z 50% -colormap jet -type uchar -omc colorslice.png

#### -oo: Output multiple images to multiple files 

Syntax: `-oo image_list` or `-oo image_spec`

Write all images on the **convert3d** stack as multiple files. There are two ways to use this command. The first is to supply a list of file names, separated by spaces: 

    c3d labelimage.nii -split -oo labelA.nii labelB.nii labelC.nii

In the above example, the image at the end of the stack will be saved as *labelC.nii*, the image next to the end of the stack will be saved as *labelB.nii* and so on. 

The second way to use the **-oo** command is to supply a pattern for the output filenames. In this case, all the images on the stack will be written. The format for the pattern is the same as for the [C++ printf command][8]. For example, the following command 

    c3d labelimage.nii -split -oo label%02d.nii

will generate images *label00.nii*, *label01.nii*, *label02.nii* and so on. The image at the end of the stack will have the highest number, and the image at the beginning of the stack will have number 00. 

#### -oomc: Output multiple multi-component images to multiple files

Syntax: `-oomc n_comp image_list` or `-oomc n_comp image_spec`

Write all images on the **convert3d** stack as multiple multi-component image files. The command is a mixture of the **-omc** and **-oo** commands. There must be a multiple of 'n_comp' images on the stack. Every consecutive 'n_comp' images on the stack will be written to a separate multi-component image.

### Commands: Stack Manipulation and Flow Control

These commands are used to manipulate the **convert3d** stack. The stack is a linear array of images. Every time an image is specified on the command line, it is loaded and placed at the end of the stack. Most operations take one image from the end of the stack, apply some operation to it, and place the result on the end of the stack. Certain commands like **-levelset** and **-reslice-matrix** take two images from the end of the stack as the input and replace them with a single output. Some other commands, like **-mean** and **-vote** take all images on the stack and replace them with a single output. 

Sometimes, for complex operations, it is useful to change the order of the images on the stack, to duplicate images, or to execute the same command multiple times. The stack manipulation and flow control commands allow you to complete complex tasks without saving intermediate images to the disk. 

#### -accum, -endaccum: Accumulate operations over all images

Syntax: `-accum command-list -endaccum`

Apply a binary operation (such as addition or multiplication) to all the images on the stack in a cumulative fashion. The command(s) will be applied to the last and second-to-last images on the stack, then to the result of this operation and the third-to-last image on the stack and so on. Below is the example of using the command to add multiple images. 

    c3d image*.nii -accum -add -endaccum -o sum.nii

#### -as: Assign image at the end of the stack to a variable

Syntax: `-as var`

Associates the image currently at the end of the stack with variable name 'var'. This allows you to retrieve the image later on the command line using the **-push** command. The **-as** and **-push** commands are useful when you need to use a certain image more than once during a convert3d operation. For example, if you want to compute the distance transform of a binary image and mask it so that the values outside of the binary image region have value 0, you would use the following command: 

    c3d binary.img -as A -sdt -push A -times -o masked_distance.img

#### -clear: Clear the image stack

Syntax: `-clear` 

Clears the image stack. Images assigned a name with the **-as** command will remain in memory. 

#### -foreach, -endfor: Loop commands over all images on the stack

Syntax: `-foreach commands-list -endfor`

This command forces the commands between **-foreach** and **-endfor** to be applied to every image on the stack. The main use of this command is to automate processing of multiple datasets. For example, 

    c3d epi*.nii -foreach -smooth 3mm -endfor -oo epism%03d.nii

#### -foreach-comp, -endfor: Loop commands over components of a multi-component image 

Syntax `-foreach-comp <N> commands-list -endfor`

This command runs the list of commands separately for each component of a set of multi-component images loaded with -mcs. This makes it possible to perform component-wise operations on multi-component images. For example, it can be used to average several multi-component images. If the image stack contains images *x1* *y1* *z1* *x2* *y2* *z2*, then the operations will be run on *[x1,x2]*, *[y1,y2]*, *[z1,z2]*. For example, if multi_1.nii to multi_10.nii are three-component images, then the mean three-component image is given by

    c2d -mcs multi_*.nii -foreach-comp -mean -endfor -omc multi_mean.nii

#### -insert: Insert image anywhere in the stack

Syntax: `-insert var pos` 

This command is similar to **-push**, but allows you to insert the image associated with 'var' at any position in the stack, counting from the end. When 'pos' is 0, the image is placed at the end of the stack (same as **-push**). When pos is one, the image is placed at the next-to-end position, and so on. 

#### -pop: Remove last image from the stack

Syntax: `-pop` 

Removes the last image from the image stack. Images assigned a name with the **-as** command will remain in memory. 

#### -popas: Remove last image from the stack and assign to variable

Syntax: `-popas var` 

Removes the last image from the stack, but also assigns it the name 'var', keeping the image in memory. Same as calling **-as** *var* followed by **-pop**. 

#### -push: Place variable at the end of the stack

Syntax: `-push var`

Places the image associated with variable name 'var' on end of the image stack. Variable names are assigned using the **-as** command. The **-as** and **-push** commands are useful when you need to use a certain image more than once during a convert3d operation. For example, if you want to compute the distance transform of a binary image and mask it so that the values outside of the binary image region have value 0, you would use the following command: 

    c3d binary.img -as A -sdt -push A -times -o masked_distance.img

#### -reorder: Rearrange images on the stack

Syntax: `-reorder k` or `-reorder fraction`

Rearranges images in the stack, such that images that are k positions apart become next to each other on the stack. In other words, if the original order of the images is 1, 2, ..., n, the new order of the images becomes 1, 1+k, 1+2k, ..., 2, 2+k, 2+2k, ..., k, 2k, ... n. Of course, n must be divisible by k. As an alternative to specifying k, you can specify a floating point number (i.e., **-reorder** 0.5), in which case k is obtained by multiplying n by the floating point number and rounding to the nearest integer. 

The following three commands are equivalent:

    c3d a1.nii a2.nii a3.nii a4.nii b1.nii b2.nii b3.nii b4.nii -reorder 4 ...
    c3d a1.nii a2.nii a3.nii a4.nii b1.nii b2.nii b3.nii b4.nii -reorder 0.5 ...
    c3d a1.nii b1.nii a2.nii b2.nii a3.nii b3.nii a4.nii b4.nii ...

The **-reorder** command us useful when you specify two sets of images using wildcards and then want to perform pairwise operations on the images. For example 

    c3d weight*.nii gray*.nii -reorder 0.5 -weighted-sum-voxelwise -o wsum.nii

is equivalent to the command

    c3d weight1.nii gray1.nii weight2.nii gray2.nii ... -weighted-sum-voxelwise -o wsum.nii

#### -dup: Duplicate the last image on the stack

Syntax: `-dup` 

Duplicates the image at the end of the stack. This is equivalent to **-as var -push var**, but shorter. An example is when you want to pass an image as both arguments to a binary operator, e.g., computing the square of the image intensity: 

    c3d input.img -dup -times -o square.img

#### -pick: Select one or more images from the stack by index

Syntax: `-pick [indices]` 

Replaces the image stack with one or more selected images. Images can be specified using 0-based indexing (0 is the first image on the stack, etc) or negative indexing (-1 is the last image on the stack). 

    c3d -mcs rgb.nii.gz -pick -1 -o blue.nii.gz
    c3d -mcs rgb.nii.gz -pick 0 -o red.nii.gz
    c3d -mcs pick -1 -2 -3 -omc bgr.nii.gz

### Commands: Voxelwise Calculations

#### -add: Voxelwise image addition

Syntax: `-add`

Adds the last two images on the stack, and places the sum at the end of the stack.

    # Add two images: x = a + b
    c3d a.img b.img -add -o x.img

    # Add three images, x = (a + b) + c in the first example, x = a + (b + c) in the second
    c3d a.img b.img -add c.img -add -o x.img
    c3d a.img b.img c.img -add -add -o x.img

    # Subtract two images, using -scale command: x = a - b
    c3d a.img b.img -scale -1 -add -o x.img

#### -atan2: Voxelwise angle from sine and cosine

Syntax: `-atan2`

Computes the angle in radians from images containing sine and cosine. This is a voxel-wise operation. It requires two images on the stack (sine followed by cosine): 

    c3d sin_theta.nii.gz cos_theta.nii.gz -atan2 -o theta.nii.gz

#### -ceil: Round up image intensities

Syntax: `-ceil `

Each image intensity is replaced by the smallest integer larger or equal to it

    c3d input.img -ceil -o output.img

#### -clip: Clip image intensity to range

Syntax: `-clip iMin iMax`

Clips image intensities, so that the values below *iMin* are converted to *iMin* and values greater than *iMax* are converted to *iMax*. This is useful for eliminating hyperintensities in images. Values *iMin* and *iMax* are intensity specifications (see below). 

    c3d mri.img -clip 1000 8000 -o mriclip01.img          // Clips below and above
    c3d mri.img -clip -inf 8000 -o mriclip02.img          // Clips above only
    c3d mri.img -clip -inf 95% -o mriclip03.img           // Clips at 95th percentile

#### -cos: Voxelwise cosine 

Syntax: `-sin`

Replaces the last image on the stack with the cosine trigonometric operation applied to all voxels. Input must be in radians.

#### -divide: Voxelwise image division    

Syntax: `-divide`

Divides one image by another. For instance to compute C = A / B, use the command 

    c3d A.img B.img -divide -o C.img

Divison may generate infinite and not-a-number (NaN) values if B contains zeros. You can use **-replace** to get rid of these values

    c3d A.img B.img -divide -replace inf 1000 -inf -1000 NaN 0 -o C2.img

#### -exp: Voxelwise natural exponent

Syntax: `-exp`

Computes exponent of each voxel in the last image on the stack.

    c3d input.img -exp -o output.img

#### -erf: Standard error function

Syntax: `-erf mu sigma`

Computes the standard error function. This is useful for applying soft thresholds. The function computes y = erf((x - mu)/sigma). 

    c3d input.img -erf 5 2 -o erf.img

#### -floor: Round down image intensities

Syntax: `-floor `

Each image intensity is replaced by the largest integer smaller or equal to it.

    c3d input.img -floor -o output.img

To round each intensity to the closest integer, use

    c3d input.img -shift 0.5 -floor

#### -log, -ln: Voxelwise natural logarithm

Syntax: `-log`

Computes natural logarithm of each voxel in the last image on the stack.

#### -log10: Voxelwise base 10 logarithm

Syntax: `-log10`

Computes base 10 logarithm of each voxel in the last image on the stack.

#### -max: Voxel-wise maximum of two images

Syntax: `-max`

Computes the voxel-wise maximum of two images. Can be used with the **-accum** command to compute maximum of all images. 

    c3d i1.nii i2.nii -max -o max12.nii
    c3d i1.nii i2.nii i3.nii i4.nii -accum -max -endaccum -o max1234.nii

#### -min: Voxel-wise minimum of two images

Syntax: `-min`

Computes the voxel-wise minimum of two images. Can be used with the **-accum** command to compute minimum of all images. 

    c3d i1.nii i2.nii -min -o min12.nii
    c3d i1.nii i2.nii i3.nii i4.nii -accum -min -endaccum -o min1234.nii

#### -mean: Mean of all images on the stack    

Syntax: `-mean `

Computes the mean of all the images on the stack. All images on the stack are replaced with the mean image.

    c3d image_*.nii -mean -o mean.nii

#### -multiply, -times: Multiply two images

Syntax: `-multiply`

Multiply two images voxel-by-voxel. The operation is applied to the last two images on the stack. 

    # Compute x = a * b
    c3d a.img b.img -multiply -o x.img

    # Compute x = a * (b + c) using add and -multiply
    c3d a.img b.img c.img -multiply -add -o x.img

Combine with the **-dup** command to compute voxelwise square of the image

    # Compute x = a^2
    c3d a.img -dup -multiply -o x.img

#### -noise-gaussian, -noise: Apply additive Gaussian noise

Syntax: `-noise-gaussian <sigma>`

Adds Gaussian noise to an image with zero mean and standard deviation sigma. Please see [Noise simulation article][15] by G. Lehmann for details.

    c3d image.nii -noise-gaussian 5 -o noisy.nii

#### -noise-poisson: Apply Poisson noise
Syntax: `-noise-poisson <scale>`

Applies Poisson (shot) noise to an image with given scale. Please see [Noise simulation article][15] by G. Lehmann for details.

    c3d image.nii -noise-poisson 5 -o noisy.nii

#### -noise-salt-pepper: Apply salt and pepper noise
Syntax: `-noise-salt-pepper <probability>`

Applies salt and pepper noise to an image with given probability. Please see [Noise simulation article][15] by G. Lehmann for details.

    c3d image.nii -noise-salt-pepper 0.1 -o noisy.nii

#### -noise-speckle: Apply Poisson noise
Syntax: `-noise-speckle <sigma>`

Applies Speckle noise to an image with given standard deviation. Please see [Noise simulation article][15] by G. Lehmann for details.

    c3d image.nii -noise-speckle 5 -o noisy.nii

#### -otsu: Otsu's thresholding
Syntax: `-otsu`

Applies the classical Otsu's binary thresholding algorithm to separate image foreground from background. Returns an image of zeros (background) and ones (foreground).

    c3d image.nii -otsu -o thresh.nii

#### -reciprocal: Image voxelwise reciprocal 

Syntax: `-reciprocal `

Computes the reciprocal of an image. For instance to compute B = 1 / A, use the command 

    c3d A.img -reciprocal -o B.img

#### -replace: Replace intensities in image

Syntax: `-replace I1 J1 I2 J2 ... `

Replace intensity I1 by J1, I2 by J2 and so on. Allowed values of intensity include **nan**, **inf** and **-inf**. 

    c3d img1.img -replace 1 128 nan 0.0 -o img2.img

#### -retain-labels: Retain labels in a label image

Syntax: `-retain-labels I1 I2 ... IN`

Assuming that the input is a multi-label segmentation image, this command keeps all labels specifed in the list and replaces the remaining labels with the background value.

    c3d seg.nii -retain-labels 2 3 4 8 -o subseg.nii

#### -rgb2hsv: Convert RGB image to HSV image

Syntax `-rgb2hsv`

Takes the last three images on the stack and treats them as red, green, and blue channels. Outputs three images corresponding to hue, saturation, value. To read color images you need the ***-msc*** command.

    c3d -mcs color.png -rgb2hsv -omc hsv.png

#### -rms: Voxelwise vector norm

Syntax: `-rms`

Computes RMS (root mean square) of all images on the stack. The command takes the square of each image on the stack, adds all the squared images and takes the square root of the result. This is very useful for statistical operations. Images must have the same size. 

    c3d img1.img img2.img img3.img img4.img -rms -o rms.img

The equivalent of this command is

    c3d img1.img img2.img img3.img img4.img -foreach -dup -times -endfor \
        -accum -add -endaccum -sqrt -o rms.img

#### -scale: Scale intensity by constant factor

Syntax: `-scale <factor>`

Multiplies the intensity of each voxel in the last image on the stack by the given factor. 

    c3d img1.img -scale 0.5 -o img2.img

#### -shift: Shift image intensity by constant

Syntax: `-shift <constant>`

Adds the given constant to every voxel.

    c3d img1.img -shift 100 -o img2.img

#### -sin: Sine

Syntax: `-sin`

Replaces the last image on the stack with the sine trigonometric operation applied to all voxels. Input must be in radians.

#### -sqrt: Take square root of image

Syntax: `-sqrt `

Computes square root of each voxel in the image.

    c3d input.img -sqrt -o output.img

#### -stretch: Stretch image intensities linearly

Syntax: `-stretch <u1 u2 v1 v2> `

Stretches the intensities in the image linearly, such that u1 maps to v1 and u2 maps to v2. The linear transformation is applied to all intensities in the image, whether inside the range or not. For example, to map a floating point image with intensities in interval (0,1) to the full range of an unsigned short image, use 

    c3d input.img -stretch 0.0 1.0 0 65535 -type ushort -o output.img

#### -thresh, -threshold: binary thresholding

syntax: `-thresh <u1 u2 vin vout> `

thresholds the image, setting voxels whose intensity is in the range [u1,u2] to vin and all other voxels to vout. values *u1* and *u2* are intensity specifications (see below). this means that you can supply values **inf** and **-inf** for u1 and u2 to construct a one-sided threshold. you can also specify *u1* and *u2* as percentiles. 
    c3d in.img -threshold -inf 128 1 0 -o out.img
    c3d in.img -threshold 64 128 1 0 -o out.img
    c3d in.img -threshold 20% 40% 1 0 -o out.img

#### -voxel-sum: Print sum of all voxel intensities

Syntax: `-voxel-sum `

Print the sum of all voxels in the image. 

    $ c3d image.img -voxel-sum 
    Voxel Sum: 200923123

#### -voxel-integral: Print volume integral of all voxel intensities

Syntax: `-voxel-integral`

Like **-voxel-sum**, but multiplies the sum of voxel intensities by voxel volume. This is useful for computing volumes of objects represented by binary images. The result is in 'ml'. 

    $ c3d image.img -voxel-integral 
    Voxel Integral: 2341

#### -wsum, -weighted-sum: Weighed sum of images with constant weights

Syntax: `-wsum weight1 weight2 ... weightN `

Computes weighted sum of the last N images on the stack. 

    c3d image1.nii image2.nii image3.nii -wsum 0.2 0.7 0.1 -o wsum.nii

This command is particularly useful for combining components in a multicomponent image. For example, for an RGB image, we can convert it to grayscale (using [ImageMagick][13] formula) as follows: 

    c3d -mcs rgb.nii -wsum 0.29900 0.58700 0.11400 -o gray.nii

#### -wsv, -weighed-sum-voxelwise: Weighed sum of images with spatially varying weights

Syntax: `-wsv `

Computes weighted sum of N weight images and N scalar images. The images must be interleaved on the stack. All images on the stack are used.

    c3d weight1.nii image1.nii weight2.nii image2.nii weight3.nii image3.nii -wsv -o mysum.nii.gz

The **-reorder** command can simplify loading the images:

    c3d weight*.nii image*.nii -reorder 0.5 -wsv -o mysum.nii.gz

### Commands: Image Header Manipulation

#### -copy-transform: Copy header information 

Syntax: `-copy-transform`

Copies the image header, specifically the image to physical space transform (origin, spacing, direction cosines), from the first image (reference) to the second image (target). This is best done with NIFTI images, which store this information well. In the example below, *out.nii* will have the same header as *first.nii* and the same intensities as *second.nii*.

    c3d first.nii second.nii -copy-transform -o out.nii

#### -mbb, -match-bounding-box: Match bounding box of one image to another

Syntax: `-mbb`

Given two images on the stack (reference and target), sets the header of the target image so that the two images occupy the same physical space. The direction cosines of the target image are set to match the reference image.  This command is related to '-copy-transform' but supports images of different size.

    c3d reference.nii target.nii -mbb -o out.nii

#### -orient: Change image orientation

Syntax: `-orient CODE`

Set the orientation of the image using one of 48 canonical orientations. The orientation describes the mapping from the voxel coordinate system (i,j,k) to the physical coordinate system (x,y,z). In the voxel coordinate system, i runs along columns of voxels, j runs along rows of voxels, and k runs along slices of voxels. It is assumed (by the NIFTI convention) that the axes of the physical coordinate system run as follows: x from (L)eft to (R)ight, y from (P)osterior to (A)nterior, z from (I)nferior to (S)uperior. 

The CODE passed in is a three-letter code consisting of letters RLAPSI. Each letter describes the anatomical direction corresponding to the voxel coordinates (i,j,k). For example, code RAI means that i runs from Right to Left, j from Anterior to Posterior, and k from Inferior to Superior. 

    c3d input.img -orient RAI -o output.img
    c3d input.img -orient SAL -o output.img

This command has the same behavior as the 'Reorient Image' menu option in ITK-SNAP. 

#### -origin: Set image origin

Syntax: `-origin vector `

Set the origin of the image. The origin is the world coordinate (in NIfTI coordinate space) of the center of the voxel (0,0,0) in the image. The origin should be specified in millimeters. 

    c3d input.img -origin 100x100x100mm -o output.img

#### -origin-voxel: Assign image origin to a voxel

Syntax: `-origin-voxel vector `

Set the origin of the image by specifying the voxel coordinates of the center of the patient (RAS) coordinate system. The vector should be specified in voxel units. 

    c3d input.img -origin-voxel 60x70x35 -o output.img
    c3d input.img -origin-voxel 50% -o output.img        # image centered around origin

#### -origin-voxel-coord: Set coordinate of specified voxel

Syntax: `-origin-voxel-coord: <index> <vector>`

This command updates the origin of the image such that the specifed voxel has the specified coordinate. For example, you can use the command to set the world coordinate (in NIFTI coordinate space) of the center voxel of the image, as follows:

    c3d input.nii -origin-voxel-coord: 50% 10x10x5mm -o output.nii

#### -set-sform: Set the transform to physical space

Syntax: `-set-sform <sform.mat> `

Sets the Nifti sform of the last image on the stack to the 4x4 matrix provided. 

#### -spacing: Set voxel spacing

Syntax: `-spacing <vector> `

Sets the voxel spacing of the image. This should always be a vector with positive components. For example, to set the spacing of the image to 1mm isotropic, use the command below. This command only changes the header of the image, not its contents. 

    c3d img.nii -spacing 1x1x1mm -o out.img

#### -swapdim: Reorder the coordinate axes of an image

Syntax `-swapdim <code>`

This command reorders the image axes (columns, rows, slices) to achieve a desired transformation between voxel space and physical space. The image remains exactly the same in physical space, but the encoding of the voxels in memory and on disk is changed to obtain the desired transformation. The transformation is specified as a three-letter 'RAI' code, as in the '''-orient''' command.

    c3d img.nii -swapdim ASL -info -o out.nii


### Commands: Image Processing

The following commands invoke an action that is applied to images. Unary commands apply the action to the last image on the stack, binary commands apply to the last two images and so on. Commands are affected by options, which are listed separately. 

#### -ad, -anisotropic-diffusion: Perona-Malik anisotropic diffusion filter

Syntax: `-ad conductance n_iter`

Executes the Perona-Malik anisotropic diffusion algorithm on the image. This smoothes the image, but with edge preservation. *Conductance* is a number between 0 and 1 that determines how well edges are preserved. *n_iter* is the number of iterations, which affects the scale of the smoothing. 

    c3d x.img -ad 0.1 100 -o ad.img

#### -biascorr: Automatic MRI bias field correction

Syntax: `-biascorr`

Performs automatic bias field correction for MRI images. This feature uses the [N3 implementation in ITK by Dr. Tustison][4], based on the N3 algorithm by Sled et al. 

    c3d mri.nii.gz -biascorr -o mricorr.nii.gz

#### -alm, -align-landmarks: Align images based on landmark matching

Syntax: `-alm dof outfile`

Performs rigid or affine alignment between to sets of landmark images. A landmark image is an image where for every intensity value, the centroid of all voxels with that intensity represents a landmark. Landmarks can be created using the paintbrush tool in ITK-SNAP (they can be spheres, cubes, etc). The first image on the stack is the target/fixed/reference image, and the second is the moving image. The parameters are the degrees of freedom, which is a number (6 for rigid, 7 for rigid+scale, 12 for affine) and the output matrix file. In this example, we have images **fixed.nii* and **moving.nii** with corresponding landmark images. We use landmarks to align the moving image to the fixed:

    c3d fixed_landmarks.nii moving_landmarks.nii -alm 6 rigid.mat
    c3d fixed.nii moving.nii -reslice-matrix rigid.mat -o moving_resliced_to_fixed.nii

#### -binarize: Convert image to binary

Syntax: `-binarize`

Converts an image to binary by mapping all background values (the background is 0 by default and can be changed by the option **-background**) to 0 and all non-background values to 1. The **-binarize** command is shorthand for the **-threshold** command. 

    c3d test.img -binarize -o binary.img 
    c3d -background 10 -binarize -o binary.img
    c3d test.img -threshold 10 10 0 1              // equivalent to above command

#### -canny: Canny edge detector

Syntax: `-canny <sigma_vector> <t_lower> <t_upper>`

Performs edge detection on the last image on the stack using the Canny filter. The parameters are a vector of standard deviations defining the scale of the edges detected and lower and upper thresholds for edge selection. See documentation on the [ITK Canny Filter][14].

#### -centroid: Report centroid of foreground voxels

Syntax: `-centroid`

Reports the centroid, in physical coordinates, of all foreground voxels in the image. 

    c3d binaryimage.img -centroid                         // centroid of all non-0 voxels
    c3d grayimage.img -thresh 1000 7000 1 0 -centroid 1   // centroid of all voxels in range 1000-7000
    c3d labelimage.img -thresh 5 5 1 0 -centroid          // centroid of all voxels with label 5
    c3d labelimage.img -split -foreach -centroid -endfor  // centroids of all labels (including 0)


#### -centroid-mark: Mark the centroid of foreground voxels

Syntax: `-centroid-mark <label>`

Marks the centroid of the foreground voxels in an image. Unlike **-centroid**, this command does not print the centroid location, but marks the closest voxel in the image with the intensity **label**. The remaining voxels are assigned 0 intensity. Combined with -dilate, this can be used to mark centers of regions with spheres.

    c3d binaryimage.nii -centroid-mark -dilate 1 3x3x3
    c3d labelimage.nii -split -foreach -centroid-mark -endfor -merge -o centers.nii


#### -color-map, -colormap: Convert scalar image to RGB using color map    

Syntax: `-color-map <ColormapName> [min max]`

Converts a scalar image to a color (RGB) image using a specified color map. The output of the command are three images, containing the red, green and blue channels of the RGB image. The mapping uses the range of the input image, e.g., using the **jet** color map, the lowest intensity pixel in the image will be mapped to blue, and the highest intesnity pixel will be mapped to red. The admissible color maps are **hot,cool,spring,summer,autumn,winter,copper,jet,hsv,red,green,blue,grey,overunder**. The command can be used with the -omc command to write RGB images. The example below generates a PNG image from a slice in a scalar image. 

    c3d scalar.nii.gz -slice z 50% -flip y -color-map jet -type uchar -omc colorslice.png

By default the full image intensity range is mapped. The optional **min** and **max** parameters can be used to set the range of the color map. 

    c3d scalar.nii.gz -slice z 50% -flip y -color-map jet 0 1 -type uchar -omc colorslice.png

#### -conv: Convolution

Syntax `-conv`

Performs convolution between the last two images on the stack. The convolution is performed using the Fourier transform. The result is an image of the same dimensions as the first image. For more details, see ["FFT Based Convolution" by Gaetan Lehmann][Lehmann].

    c3d image.nii kernel.nii -conv -o result.nii


[Lehmann]: https://hdl.handle.net/10380/3154

#### -comp, -connected-components: Compute connected components

Syntax: `-comp`

Computes the connected components of a binary image. Each connected component is assigned an integer index. Indices are ordered by the size of the component, so the component assigned index 1 is the largest. The background is assigned index 0. To select the largest connected component, combine the call to **-comp** with a call to **-threshold**. 

    c3d binary.img -comp -o comp.img
    c3d binary.img -comp -threshold 1 1 1 0 -o largest_comp.img

#### -cmv, -coordinate-map-voxel: Generate voxel coordinate maps (voxel units)

Syntax: `-cmv`

For a *N*-dimensional image, replaces the last image on the stack with *N* images. The *k*-th output image at each voxel contains the $k$-th coordinate of that voxel, in voxel units.

    c3d image.nii -cmv -oo coordmap%d.nii.gz

One can use this command to split a brain segmentation image into a left hemisphere segmentation and a right hemisphere segmentation (assuming the X coordinate corresponds to the right-left axis)

    c3d seg.nii -as SEG -cmv -pop -pop  -thresh 50% inf 1 0 -as MASK \
        -push SEG -times -o seg_left.nii.gz \
        -push MASK -replace 1 0 0 1 \
        -push SEG -times -o seg_right.nii.gz

#### -cmp, -coordinate-map-physical: Generate voxel coordinate maps (voxel units)

Syntax: `-cmp`

This command is similar to **-cmv** (**-coordinate-map-voxel**), but the output will contain the physical coordinates of the voxels, in the NIFTI (RAS) coordinate frame. 

#### -create: Generate blank image

Syntax: `-create dimensions voxel_size`

Creates a new blank image with specified dimensions and voxel size, and places it at the end of the stack. The image is set to the current background value, which is 0 by default but can be overwritten with the **-background** command. The origin of the image can be changed with the **-origin** command. 

    c3d -create 256x256x160 1x1x1mm -o newimage.img
    c3d -background 128 -create 256x256x160 1x1x1mm -origin 128x128x80mm -o newimage.img

#### -dilate: Binary dilation

Syntax: `-dilate <label> <radius_vector>`

Applies the dilation [mathematical morphology][5] operation to a binary image. The first parameter is the intensity value of the object that is to be dilated. The second is the radius of the dilation structuring element in 3D. 

    c3d binary.img -dilate 255 3x3x3vox -o newimage.img

#### -erode: Binary erosion

Syntax: `-erode <label> <radius_vector>`

Applies erosion [mathematical morphology][5] operation to a binary image. The first parameter is the intensity value of the object that is to be eroded. The second is the radius of the erosion structuring element in 3D. 

    c3d binary.img -erode 255 3x3x3vox -o newimage.img

#### -export-patches, -xp: Fixed size patch sampling from masked regions

Syntax: `-export-patches <outfile> <radius_vector> <frequency>`

See also: **-export-patches-aug (-xpa)** command, which provides data augmentation for deep learning.

This command samples patches from a region of a ND image and stores them into a data file that can be read easily in other software, for example, NumPy. This is useful for generating training data for machine learning projects. Multiple "channels" can be sampled.

    c3d chan1.nii chan2.nii chan3.nii mask.nii -xp samples.dat 4x4x4 100

This command will sample the three images chan1, chan2, chan3 at foreground voxels in the mask. Voxels in the mask foreground region are sampled randomly, following a uniform distribution. The value of 100 means that every 100-th voxel, on average, is sampled. The radius 4x4x4 means that patches of size 9x9x9 will be generated. For each sampled voxel, the sampled intensity data is represented as a 3x9x9x9 array in this example.

To read these samples in NumPy use the following code

    dims = (9,9,9)                          # Patch dimensions
    k = 3                                   # Number of channels
    bps = (4 * k * reduce(mul, dims, 1))    # Bytes per sample
    np = os.path.getsize(fname) // bps      # Number of samples
    arr = numpy.memmap(fname,'float32','r',shape=(np,k) + dims)

It is also possible to visualize the extracted samples in ITK-SNAP by reading them as a raw image, with dimensions equal to the dimensions of the patch, and the z-dimension multiplied by the number of samples.

The command can also be used to extract entire structures. For example, if we have a binary segmentation of a lesion of an approximately known size in an MRI scan, we can extract a patch of given size centered on this lesion, as follows:

    c3d mri.nii lesion_seg.nii -centroid-mark 1 -xp single_sample.dat 50x50x20 1

In the above example, **-centroid-mark** transforms the lesion segmentation into a single-voxel mask, from which the sample from the MRI is taken.

#### -export-patches-aug, -xpa: data augmentation for deep learning

Syntax `-export-patches-aug <N> <sigma_angle>`

This command must precede the `-export-patches (-xp)` command and instructs this command to not only sample patches along the image axes but to also sample **N** randomly rotated patches. Rotation is around a uniformly distributed axis with a rotation angle distributed normally with teh standard deviation **sigma_angle**, specified in degrees. This kind of sampling is useful for data augmentation for machine learning algorithms.

    c3d chan1.nii chan2.nii chan3.nii mask.nii -xpa 5 10 -xp samples.dat 4x4x4 100

#### -fft: Fast Fourier transform

Syntax `-fft`

Computes the Fourier transform of a real-valued image at the end of the stack. The image is replaced by the real and imaginary components of the FFT. This command is only available if **convert3d** is compiled with the FFTW library support. 

    c3d image.nii -fft -oo real.nii imag.nii

#### -flip: Flip image around an axis    

Syntax: `-flip axes`

Flips the image around specified axes. The parameter 'axes' may be any combination of characters 'x', 'y', and 'z'; the order does not matter. 

    c3d input.img -flip xy -o output.img

#### -glm: General linear model    

Syntax: `-glm design_matrix_file contrast_vector_file`

Applies voxel-wise general linear model to a set of images. More precisely, the general linear model solves the following system: $Y = X \beta + \epsilon$, where Y are the observations (a list of n images, where each voxel is treated as an independent observation); X is the $n x k$ design matrix, where $k$ is the number of factors; $\beta$ is a set of $k$ unknown images (factors) and $\epsilon$ is the error term. The command will compute the $\beta$ images and return a weighted sum of them, where the weights are specified in the contrast vector. The design matrix and the contrast vector are passed in as files. The file format is just a space-separated list of numbers. For a good explanation of the general linear model, see [S. Kiebel and A. Holmes, General Linear Model, in Ashburner, Friston, Holmes eds., *Human Brain Function, 2nd Edition*][6]. The example below computes the regression coefficient between a set of longitudinal images and subject's age: 

      echo "1 67.00" > design_mat.txt
      echo "1 75.00" >> design_mat.txt
      echo "1 80.00" >> design_mat.txt
      echo "1 83.00" >> design_mat.txt
      echo "0 1" >> contrast_vec.txt
      c3d time1.img time2.img time3.img time4.img -glm design_mat.txt contrast_vec.txt -o regress.img

#### -grad, -gradient: Image gradient

Syntax `-grad`

Computes the gradient of the image. Each component of the gradient is placed on the stack in order (x,y,z). The gradient is computed in physical RAS coordinates, taking into account image spacing and orientation. In other words, the gradient is the vector in physical space orthogonal to the isocontours of the image. No smoothing is performed, so it is a good idea to smooth the image first before computing the gradient.

    c3d myimage.nii -smooth 1.2vox -grad -oo grad_comp_%02d.nii

#### -hesseig, -hessian-eigenvalues: Compute eigenvalues of the Hessian matrix

Syntax `-hesseig <scale>`

Computes the Hessian matrix at every pixel of an image and the eigenvalues of the Hessian. Images of the eigenvalues (sorted by value) are placed on the stack. These images are useful as texture features. See also the '''-steig''' command. The scale determines the amount of Gaussian smoothing applied for computing the partial derivatives in the Hessian, and is in physical (mm) units.

    c3d myimage.nii -hesseig 2.0 -oo eig%02d.nii.gz

#### -hessobj, -hessian-objectness: Hessian objectness filter

Syntax: `-hessobj <dimension> <min_scale> <max_scale>`

Also known as the Frangi vesselness filter, this filter can be used to highlight tube-like, sheet-like and blob-like objects in the image. For details, see documentation to the [corresponging ITK class][HTOMIF]. 

Parameter `dimension` is an integer that determines the kind of features that are highlighed. Use 0 for blobs, 1 for tubes, 2 for pancakes, etc. The min and max scale parameters are floating point values, giving the scale of the features highlighted, in physical units. Typically, just one scale is used.

    # Detect vessel-like structures at scale 0.5mm
    c3d image.nii.gz -hessobj 1 0.5 0.5

 [HTOMIF] http://www.itk.org/Doxygen/html/classitk_1_1HessianToObjectnessMeasureImageFilter.html

#### -holefill: Fill holes in binary image

Syntax: `-holefill intensity_value [0|1] `

Apply the binary hole filling algorithm to a particular intensity value in the image. The input image is typically a binary image or a multi-label segmentation image. Holes (voxels not matching the specified intensity value that are completely contained by voxels matching this value) are filled. The second parameter specifies what type of topological connectivity is used to determine holes. The value 0 uses the default algorithm in ITK (face connectivity) and 1 uses the full connectivity variant (face, edge and vertex connectivity). For more details see the [ITK page for this algorithm][7]. 

    c3d segmentation.nii.gz -holefill 5 0 -type uchar -o filledlabel5.nii.gz 

#### -lstat, -label-statistics: Display segmentation volumes and intensity statistics

Syntax: `-lstat`

Given a grayscale image and a multilabel (or binary) image, this command computes the statistics for every label in the latter, including volumes, average grayscale intensity, etc. For instance, if image *mri.nii* is a medical image and *seg.nii* is a multilabel segmentation of the image with labels 0, 1 and 4, the following command can be used to print the statistics of the intensity of *mri.nii* for each of the labels 

    c3d mri.nii seg.nii -lstat

The output contains the mean, standard deviation, maximum intensity and minimum intensity for each label. If you just need volumes from a multi-label image, use **-dup** command as follows:

    c3d seg.nii -dup -lstat

#### -lms, -landmarks-to-spheres: Generate image of spheres from landmark data

Syntax: `-lm <landmark_file> <radius>`

Reads a list of landmarks from a file and for each landmark draws a sphere with the landmark at the center. Spheres are drawn over the last image on the stack. Landmarks are specified in a file with four entries per line: *x y z value*, where *(x,y,z)* give the landmark position in physical units (mm) and value is the image value that the sphere will be filled with. The radius (for now) is specified in physical units (mm). The example below first creates a landmark file, and then calls c3d to create the spheres. Notice that the **-scale** command is first used to create a blank image of the same dimensions as the input image. 

    echo 14.5 23.12 44.53 1 > landmarks.txt
    echo 76.3 43.12 34.32 2 >> landmarks.txt
    ...
    c3d image.img -scale 0 -lms landmarks.txt 3 -o landspheres.img

The command also supports an alternative format for the landmarks that allows specification in voxel
units and percent

    echo 50%x50%x50% 1 > landmarks.txt
    echo 25x25x30vox 1 >> landmarks.txt
    echo 10x10x45mm 2 >> landmarks.txt 
    ...
    c3d image.img -scale 0 -lms landmarks.txt 3 -o landspheres.img

#### -laplacian, -laplace: Laplacian filter

Syntax: `-laplacian`

Applies the Laplacian filter to the image. Used to detect ridges of intensity. Typically, used with the **-smooth** option to obtain the equivalent of convolving the image with the *Laplacian of the Gaussian (LoG)* kernel: 

    c3d input.img -smooth 1.2vox -laplacian -o output.img

#### -levelset: Level set segmentation

Syntax: `-levelset n_iter `

Perform level set segmentation for *n\_iter* iterations, like in ITK-SNAP. The last image on the stack is treated as the initialization image and the next-to-last image on the stack is the speed image. Both images should be in the range between -1 and 1. Here is how the signs of the different images are interpreted 

|    | Speed Image   | Initialization Image | Output Image |
| -- | ------------- | -------------------- | ------------ |
| +1 | Foreground    | Outside              | Outside      |
| -1 | Background    | Inside               | Inside       | 

Here is an example where you have the speed and the initialization given: 

    c3d speed.img initial.img -levelset-curvature 0.5 -levelset 100 -o seg.img

Here is an example of segmenting the ventricles in an MRI image, where the ventricles and other CSF have intensity below 715. The image seg_bubbles.nii.gz in this example is a binary image of the initialization seeds (1 inside the seeds, 0 outside). 

    c3d brain.nii.gz -erf 715 100 -scale -1 seg_bubbles.nii.gz \
        -replace 0 1 1 -1 -levelset-curvature 0.2 -levelset 500 \
        -thresh -inf 0 1 0 -o segmentation.nii.gz

Another example of smoothing a binary image that is useful for cleaning up manual segmentations. Here the speed image is positive inside the binary object, and the initialization is negative inside the object. The command writes out both the level set image (whose 0-level set is the smoothed boundary of the binary object) and the smoothed binary object 

    c3d binary.img -threshold 1 inf 1 -1 -binary.img 1 inf 1 -1 \
        -levelset-curvature 1.5 -levelset 100 -o levelset.img \
        -thresh -inf 0 1 0 -o smoothed_binary.img

#### -mf, -mean-filter: Mean filter

Syntax: `-mf <radius_vector>`

Applies the mean filter: the intensity of each voxel is replaced by the mean of the intensities in the neighborhood of size specified by the radius parameter. For example, the following code will apply the mean filter with the 5x5x5 neighborhood. 

    c3d in.nii -mf 2x2x2 -o filtered.nii

#### -median, -median-filter: Median filter

Syntax: `-median <radius_vector>`

Applies the median filter: the intensity of each voxel is replaced by the median of the intensities in the neighborhood of size specified by the radius parameter. For example, the following code will apply the median filter with the 5x5x5 neighborhood. 

    c3d in.nii -median 2x2x2 -o median.nii

#### -merge: Merge images from previous split command   

Syntax: `-merge`

Works in conjunction with the **-split** command. Has similar behavior to **-vote**, except that label values are carried from the input to the **-split** command. 

#### -msq, -mean-square: Compute mean square difference metric

Syntax: `-msq [movtransform.mat] [reftransform.mat]`

Compute the mean square difference metric between the last two images on the stack. If an optional *movtransform.mat* file is provided, the metric is computed by applying the transform to the moving image. If, in addition to *movtransform.mat*, the optional *reftransform.mat* file is also provided -- the moving transform is applied to the moving image, the ref transform is applied to the reference image, and the metric is computed in an image space that is physically halfway between the reference and moving images. This may be useful for unbiased metric computation if the two transforms are inverse of each other as both images undergo similar amount of interpolation. The definitions of reference and moving images and the transform file format are similar to the **-reslice-matrix** command. 

    # Compute metric between ref.nii and mov.nii
    c3d ref.nii mov.nii -msq

    # Compute metric between ref.nii and mov.nii after applying transform to mov.nii
    c3d ref.nii mov.nii -msq tmov.mat

    # Compute metric between ref.nii and mov.nii in a neutral space after applying transforms to both
    c3d ref.nii mov.nii -msq tmov.mat tref.mat

#### -mi, -mutual-info: Compute mutual informaiton metric

Syntax: `-mi [movtransform.mat] [reftransform.mat]`

Compute the mutual information metric between the last two images on the stack. See documentation for **-msq**.

#### -mmi, -mattes-mutual-info: Compute mutual informaiton metric

Syntax: `-nmi [movtransform.mat] [reftransform.mat]`

Compute the Mattes mutual information metric between the last two images on the stack. See documentation for **-msq**.

#### -ncc, -normalized-cross-correlation: Compute normalized cross-correlation image

Syntax: `-ncc <radius_vector>`

Computes normalized cross-correlation between two images that occupy the same physical space. Each voxel in the resulting image is the cross-correlation of patches of given radius surrounding the voxel in the two input images. This is different from **-ncor**, which computes a global cross-correlation metric value. 

#### -nmi, -normalized-mutual-info: Compute mutual informaiton metric

Syntax: `-nmi [movtransform.mat] [reftransform.mat]`

Compute the normalized mutual information metric between the last two images on the stack. See documentation for **-msq**.

#### -ncor, -normalized-correlation: Compute normalized correlation metric

Syntax: `-ncor [movtransform.mat] [reftransform.mat]`

    :   Compute the normalized correlation metric between the last two images on the stack. See documentation for **-msq***.

#### -nlw, -normalize-local-window: Standardize image intensity using local neighborhood

Syntax: `-nlw <radius>`

This command takes as inputs an image and a mask image. At each voxel, the mean of the local neighborhood is subtracted, and the result is divided by the standard deviation of the neighborhood. The mean and standard deviation are computed only over the masked region. You might also want to multiply by the mask.

    c3d gray.nii.gz mask.nii.gz -nlw 10x10x10 -o residual.nii.gz

#### -oli, -overlay-label-image: Overlay segmentation image on grayscale image

Syntax: `-oli lookup_table_file opacity`

This command takes a grayscale image and a label image (i.e. image with a set of discrete values) and produces red, green and blue components of a color image. The resulting color image is an overlay of the labels over the grey image. The first parameter (*lookup\_table*) is a text file with entries in the format 

    label_value red green blue alpha 

Alpha values must be between 0 and 1. Red, green and blue values should be on the same order as the intensity of the grey image (typically 0-255). The text file is compatible with ITK-SNAP and can be generated using the ITK-SNAP `Segmentation->Save Label Descriptions` command. The second parameter (*opacity*) is between 0 and 1 and sets the overall opacity of the overlay. The output of this command is similar to the way ITK-SNAP presents segmentation data on top of grayscale images. 

    c3d gray.nii.gz -stretch 2% 98% 0 255 -clip 0 255 seg.nii.gz -oli labels.txt 0.5 -omc rgb.nii.gz

Note: this command does not interpolate between entries in the lookup table. It should not be used for images with a continuous intensity spectrum. 

Here is a more complex example, used to visualize a segmentation result. We do a few things in this command: trim grayscale and segmentation images to an ROI around the object of interest; map intensity range of the grayscale image to 0-255; extract slices through the middle of the cropped images; overlay segmentation on the grayscale image; and save as a color PNG file. 

    c3d seg.nii.gz -trim 20x20x0vox -as S gray.nii.gz -stretch 2% 98% 0 255 -clip 0 255 \\
        -reslice-identity -push S -foreach -slice z 50% -flip xy -endfor \\
        -oli labels.txt 0.5 -type uchar -omc ovl.png

#### -overlap: Compute relative overlap between binary images    

Syntax: `-overlap Z`

Compute relative overlap between labels in the last two images on the stack. Overlap is computed for a given label **Z**, i.e., the number of voxels that are equal to **Z** in both images is computed and divided by either the average number of voxels equal to **Z** in both images (to get Dice coefficient) or by the size of the region where at least one of the images is equal to **Z** (Jaccard coefficient). 

The command below computes overlap for label 255.

    c3d -verbose seg1.img seg2.img -overlap 255

The output of the command is in the following terse format, with the last two values giving Dice and Jaccard coefficients, respectively. 

    OVL: 1, 2383, 2474, 1807, 0.744081, 0.592459

Use the flag **-verbose** to get full information.

    Matching voxels in first image:  2383
    Matching voxels in second image: 2474
    Size of overlap region:          1807
    Dice similarity coefficient:     0.744081
    Intersection / ratio:            0.592459

This command does not alter the stack.

#### -pad: Pad image with constant value

Syntax: `-pad <padlower> <padupper> <value> `

Pads the image by a given percentage or number of voxels. The *padlower* dimension pads along the zero faces of the image, and the *padupper* dimension pads along the upper faces of the image. For example to add 1 voxel to the left side of an image, do 

    c3d img1.nii -pad 1x0x0vox 0x0x0vox 0 -o padded.nii

while 

    c3d img1.nii -pad 2x2x4vox 0% 0 -o padded.nii

adds two voxels padding to the left and posterior sides, and four slices to the bottom of the image. Note that the first argument changes the location of voxel (0,0,0) and thus the origin of the output image will be changed to maintain anatomical alignment between the padded and original images. 

Normally you will want to pad with zeros, but you can pad with any constant value, eg : 

    c3d img1.nii -pad 10% 10% 1 -o padded.nii

Adds 10% to all sides of the image, and fills the new voxels with the value 1. 

#### -pad-to: Pad image to desired size with constant value

Syntax: `-pad <target_size> <value> `

Pads the image symmetrically with constant value to achieve a desired size.

    c3d image1.nii -pad-to 60x60x100 0 -o padded.nii


#### -pca: Principal components analysis of foreground voxels

Syntax: `-pca`

Similar to the *-centroid* command, computes the centroid and prinicipal components of the foregrond voxels in the image. For example if the image is a binary image of an ellipsoid, this will report the center and the principal axes of the ellipsoid, in physical NIFTI coordinates.

    c3d binaryimage.img -pca                              // centroid of all non-0 voxels
    c3d labelimage.img -thresh 5 5 1 0 -pca               // centroid of all voxels with label 5
    c3d labelimage.img -split -foreach -pca -endfor       // centroids of all labels (including 0)

#### -probe: Report image intensity at a voxel

Syntax: `-probe <point_spec>`

Prints the value of the image at the position specified by the parameter `point_spec`, which may be in physical units or voxel units:

    c3d img1.img -probe 128x120x160vox
    c3d img1.img -interpolation NearestNeighbor -probe 60x60x60mm
    c3d img1.img -probe 50%

#### -rank: Voxelwise ranking of intensity values

Syntax: `-rank `

This command takes N images as the input (all the images on the stack are used). It also generates N images as the output. For voxel k in image j, it assigns it a label based on its rank among the values of voxel k in all N images. If the voxel has highest intensity in image j, then the j'th output will have value 1. 

    c3d img1.img img2.img ... imgN.img -rank -oo rank%d.img

#### -rf-train: Train Random Forest classifier

Syntax: `-rf-train <classifier_file>`

This command trains a classifier using an implementation of the [Breyman et al. Random Forest Algorithm][Br2001], with modifications proposed by [Criminisi and Shotton][Cr2004]. The stack must contain one or more images of features (e.g., grayscale images), followed by a multi-label image. The latter must have at least two non-zero labels corresponding to different classes. The classifier is trained on a voxel by voxel basis. All voxels with label *L* are treated as the examples of class *L*. The classifier is output to a binary file that can later be used by the **-rf-apply** command. Multiple parameters can be specified with the **-rf-param-xxx** options before calling **-rf-train**. The stack is not modified by this command.

    # Training with two MRI modalities as features and default parameters
    c3d t1_mri.nii t2_mri.nii segmentation.nii -rf-train myforest1.rf

    # Training with patches as features (see docs for -rf-param-patch)
    c3d ultrasound.nii seg.nii -rf-param-patch 2x2x2 -rf-train myforest2.rf

    # Applying the classifier
    c3d ultrasound.nii -rf-apply myforest2.rf -omc class_prob.nii.gz

The commands are meant to replicate the "classification" pre-segmentation mode in ITK-SNAP, i.e., extending a rough example segmentation to the entire image domain. It is possible to also use the commands to train classifiers jointly on data from multiple subjects, each with its own segmentation, as long as the images from the different subjects occupy the same image space and can be stacked into a 4-dimensional image. For example:

    # Train using MRI and segmentations from N subjects
    c4d mri_subj*.nii -tile w -popas ALLMRI \
        seg_subj*.nii -tile w -popas ALLSEG \
        -rf-param-patch 2x2x2x0 \
        -push ALLMRI -push ALLSEG -rf-train myforest.rf

    # Apply using single MRI
    c4d mri_new.nii -rf-apply myforest.rf -omc classprob.nii

 [Br2001] Breiman, L. (2001). Random forests. Machine learning, 45(1), 5-32.
 [Cr2004] Criminisi, A., & Shotton, J. (2013). Decision forests for computer vision and medical image analysis. Springer Science & Business Media

#### -rf-apply: Apply Random Forest classifier 

Syntax: `-rf-apply <classifier_file>`

This command applies a classifier trained previously by **-rf-train**. The stack must contain the same number of feature images as when training. The images will be removed from the stack and replaced with a set of K probability images, where K is the number of classes during training. See examples under **-rf-train** for usage.

#### -region: Extract region from image

Syntax: `-region vOrigin vSize `

Extract a rectangular region from the image. The first parameter is the position of the corner of the region, and the second is the size of the region. 

    c3d img1.img -region 20x20x20vox 50x60x70vox -o img2.img
    c3d img1.img -region 25% 50% -o img3.img

#### -resample: Resample image to new dimensions

Syntax: `-resample <dimensions> `

Resamples the image, keeping the bounding box the same, but changing the number of voxels in the image. The dimensions can be specified as a percentage, for example to double the number of voxels in each direction. The **-interpolation** flag affects how sampling is performed. 

    c3d img1.img -resample 123x142x200 -o img2.img 
    c3d img1.img -resample 200% -o img2.img 
    c3d img1.img -resample 100x100x200% -o img2.img 
    c3d img1.img -background 4.0 -interpolation Cubic -resample 123x142x200 -o img2.img

#### -resample-iso: Resample image to (approximately) isotropic resolution

Syntax: `-resample-iso <min|max>`

Resamples the image to have approximately isotropic resolution, either based on the smallest voxel dimension ('min' mode) or largest voxel dimension ('max' mode). This command calls **-resample** with appropriately calculated new image dimensions. The bounding box of the image in physical space is preserved. Therefore, since the image dimensions must be integer, the actual voxel dimensions after resampling may not be precisely isotropic. 

    c3d img1.img -resample-iso min -o img2.img


#### -resample-mm: Resample image to new resolution

Syntax: `-resample-mm <voxel_size> `

Resamples the image as in **-resample**, but the user specifies the new voxel size rather than dimensions. This may not be precise, so the bounding box of the image may change. A warning will be generated in that case. 

    c3d img1.img -resample-mm 1.0x1.5x1.5mm -o img2.img

#### -reslice-matrix: Resample image using affine transform

Syntax: `reslice-matrix <transform_file>`

Applies affine transform to an image. The first image on the command line is the reference image. The output image will have same dimensions and voxel-to-space transform as the reference image. The second parameter is the moving image, i.e., the image to be resliced. The transform file is a 4x4 matrix that gives a transform *from the physical space of the reference image to the physical space of the moving image*. To generate a transform file, you must use some sort of registration algorithm, which is not part of **convert3d** (at least not yet).

    c3d reference_image.nii moving_image.nii -reslice-itk transform_file.txt -o output_image.nii

There is a separate tool **c3d_affine_tool** for manipulating these matrix files.

Here is an example of using **-reslice-matrix** with SPM. In MATLAB, you would execute the following commands 

    p = spm_coreg('reference.nii', 'moving.nii');
    m = inv(spm_matrix(p(:)'));
    save -ascii ref2mov_mat.txt

Then to reslice using **convert3d**, run 

    c3d reference.nii moving.nii -reslice-matrix ref2mov.txt -o resliced.nii

The related command **-reslice-itk** can be used to use affine transforms computed by ANTs software.

#### -reslice-itk: Resample image using affine transform

Syntax: `-reslice-itk <transform_file> `

Applies affine (or other) transform in ITK (ANTs) format to an image. See notes to **-reslice-matrix** for usage.

#### -reslice-identity: Resample image using identity transform 

Syntax: `-reslice-identity `

Applies the **-reslice-matrix** command with the identity transform. This is useful when you have two scans of the same subject with different coordinate transformations to patient space and you want to resample one scan in the space of another scan. For example, if you have T1 and T2 images in different coordinate frames, and want to reslice the T2 image into the space of the T1 

    c3d t1.nii t2.nii -reslice-identity -o t2_in_t1_space.nii

#### -sdt, -signed-distance-transform: Signed distance transform of a binary image

Syntax: `-sdt`

Computes the signed distance transform of a binary image. Voxels where the binary image is non-zero will have negative values and voxels where the binary image is zero will have negative values. The magnitude of the value will be the approximate Euclidean distance to the boundary of the object represented by the binary image. 

    c3d binary.img -sdt -o dist.img 

#### -sharpen: Sharpen edges in the image

Syntax: `-sharpen`

Applies the Laplacian sharpening filter from ITK, which accentuates the edges in the image.

    c3d input.nii.gz -sharpen -o output.nii.gz

#### -slice: Extract slices from an image

Syntax: `-slice axis position_spec`

Extracts a slice along the specified axis (x,y or z). The position specifier **position_spec** can be a single slice or a range of slices. For a single slice, it can be specified as a number or a percentage. Numbering is zero-based, i.e, the first slice is slice 0, the last slice is N-1, where N is the number of slices. For a range, use MATLAB notation first:step:last. The slice is placed on the stack as an image with size 1 in the last dimension. You can save the slice as a 2D PNG image. 

    c3d input.img -slice x 128 -o myslice.nii.gz
    c3d input.img -slice y 50% myslice.nii.gz
    c3d input.img -slice z 25% -type uchar -stretch 0 2000 0 255 -o myslice.png
    c3d input.img -slice z 0:-1 -oo slice%0d.nii.gz 
    c3d input.img -slice z 20%:10%:80% -oo slice%0d.nii.gz 

With the new command **c4d**, the **-slice** command can be used to extract volumes from a 4D image. This can be useful to reformat a 4D NIFTI image as a 3D multi-component NIFTI image, using the command 

    c4d input4d.nii.gz -slice w 0:-1 -omc output3d_multicomp.nii.gz

#### -slice-all: Extract slices from all images on the stack

Syntax `-slice-all axis position_spec`

This command behaves identical to the **-slice** command, but all images on the stack are sliced, and the slices are interleaved. This is useful for slicing multi-component images. For example, if you read a four-component image 'test4.nii.gz', you can extract and save the slices as follows:

    c3d test4.nii.gz -slice-all 20%:10%:80% -oomc 4 slice4_%03d.nii.gz

#### -smooth: Gaussian smoothing

Syntax: `-smooth <sigma_vector> `

Applies Gaussian smoothing to the image. The parameter vector specifies the standard deviation of the Gaussian kernel. Also see [Vector Format Specification][10] below. 

    c3d img1.img -smooth 2x1x1vox -o out.img

#### -smooth-fast: Fast approximate Gaussian smoothing

Syntax: `-smooth-fast <sigma_vector> `

Applies Gaussian smoothing to the image using the fast [Deriche recursive smoothing algorithm][15].  The parameter vector specifies the standard deviation of the Gaussian kernel. Also see [Vector Format Specification][10] below. 

    c3d img1.img -smooth-fast 20x10x10vox -o out.img

#### -split: Split multi-label image into binary images

Syntax: `-split`

This command takes a multilabel image (one with a small number of discrete intensity levels), and replaces it with a set of binary images, one for each of the levels. The images can later be recombined using the **-merge** command. The labels corresponding to each binary image are remembered by **convert3d** so that when **-merge** is called, the labels are faithfully reassigned. The **-merge** command treats each input as a probability image, and selects at each voxel the label that has highest probability. The example below smooths each label independently, then recombines using **-merge** 

    c3d multilabel.nii -split -foreach -smooth 3mm -endfor -merge -o ml_smooth.nii

Also of note is that the **-split** command will disregard infinite intensity values. So if you want to apply voting to a subset of the labels, you can replace labels you do not care about with *inf*, for example, using the **-thresh** command. 

#### -staple: STAPLE algorithm to combine segmentations

Syntax: `-staple <intensity_value> `

Runs the ITK implementation of the STAPLE algorithm ([See Paper][11]). STAPLE generates an estimate of the 'true' segmentation of a structure given a set of segmentations by different raters. This command treats all images on the stack as inputs. Each image is considered to be a segmentation by a different rater. The parameter *intensity_value* specifies the label in the segmentation images corresponding to the structure of interest (e.g., the segmentation image may have value 1 corresponding to the caudate and value 2 corresponding to the hippocampus. To run STAPLE on the hippocampus, pass in 2 as the *intensity_value*). The output of STAPLE is a real-valued image with voxels between 0 and 1, representing the probability of each voxel being in the 'true' segmentation. This image can be thresholded to get a binary consensus segmentation. Additional outputs (estimates of the sensitivity and specificity of each rater) are printed out if the **-verbose** command is used before the **-staple** command. 

    c3d -verbose rater1.img rater2.img rater3.img -staple 1 -o probmap.img
    c3d -verbose rater*.img -staple 1 -threshold 0.5 inf 1 0 -o bin_segm.img

#### -steig, -structure-tensor-eigenvalues: Compute eigenvalues of the structure tensor

Syntax `-steig <scale> <radius>`

Computes the Hessian matrix at every pixel of an image and the eigenvalues of the Hessian. Images of the eigenvalues (sorted by value) are placed on the stack. These images are useful as texture features. See also the '''-steig''' command. The scale determines the amount of Gaussian smoothing applied for computing the partial derivatives in the Hessian, and is in physical (mm) units.

    c3d myimage.nii -hesseig 2.0 -oo eig%02d.nii.gz

#### -test-image, -test-probe: Test condition

Syntax: `-test-image [tolerance]` and `-test-probe <vector> <value> [tolerance]`

These advanced commands (with more to come in the future) are primarily meant to allow testing of **c3d**. However, they can also be used for flow control in shell scripts (e.g., **bash** shell). The commands check a certain aspect of the **c3d** state and cause the program to exit with either return code 0 if the test succeeded or a non-zero return code if the test failed. 

**-test-image** tests if the last two images on the stack are identical (both in terms of data and header). Returns 0 if the images are identical. The optional tolerance parameter has default value 1e-8. 

    c3d input1.img input2.img -test-image

**-test-probe** is similar to the **-probe** command. It tests if the value of the last image on the stack at the position given by **vector** is equal to the **test_value**. An optional tolerance value may be specified, the default is 1e-8. 

    c3d input1.img -test-probe 40x40x20vox 1.0 1e-6

#### -tile: Tile and stack multiple images into one

Syntax: `-tile <tile_spec>`

Tiles multiple images into a single image -- including stacking slices into a 3D volume. The command takes all images on the stack and produces a single tiled image. The **tile_spec** parameter can either specify a coordinate axis (x, y, or z) along which to tile the images, or a layout vector (e.g., **4x4**) which specifies the tiling along each coordinate. Passing 0 for the last value in the layout vector determines the value based on the number of images currently loaded. For example, to create a 3D volume from a set of slices, we use 

    c3d slices*.png -tile z -o volume.nii.gz

And to arrange the same 2D slices into a 2D montage of 4 images per row, we would use the **c2d** command as follows: 

    c2d slices*.png -tile 4x0 -type uchar -o montage.png

#### -trim: Trim background region of image

Syntax: `-trim <margin_vector>`

Use this command to trim background in an image. When most of the image is filled by background, this command will find the smallest rectangular region that contains all of the non-background voxels in the image. I will then expand this region by the margin of the size specified, and return the resulting region as the new image. For example, this command will trim an image, leaving a 5-voxel margin of background values on all sides

    c3d in.img -trim 5vox -o out.img

#### -trim-to-size: Trim image to given size

Syntax: `-trim-to-size <size_vector>`

Like **-trim**, this command trims the background in an image. However, instead of **-trim**, you specify the target size of the output region. The actual region may be smaller if the specified region falls outside the boundaries of the input image. For example, if you want a 64x64x128 image containing all the foreground pixels in your image, call 

    c3d in.img -trim-to-size 64x64x128vox -o out.img

#### -vote: Vote among images on the stack

Syntax: `-vote `

This command takes all images on the stack as arguments and at each voxel *(i,j,k)* returns the index of the image for which the image value at *(i,j,k)* is the greatest. This is most useful when combining probability maps into a single label image. If images prob1.img, prob2.img, etc. give the probability of label 1, 2, etc. over the image domain, the **-vote** command will return the most probable label at each voxel. 

    c3d prob1.img prob2.img prob3.img -vote -type uchar -o label.img

The value assigned to each image is based on its position from the bottom of the stack, with zero indicating bottom-most image. In the example above, the output image has values 0 for voxels where prob1.img is highest, 1 for prob2.img and 2 for prob3.img. Also see the related commands **-split** and **-merge**. 

#### -vote-mrf: Vote with Markov Random Field regularlization

Syntax: `-vote-mrf <mode> <lambda>`

This command is similar to **-vote** but it performs regularlization using the Markov Random Field (MRF). This form of regularization penalizes the total surface area of the segments in the output. It results in more contiguous segments. 

The command takes all the images on the stack and assumes that they are likelihood images corresponding to labels 1, 2, ... N. This means that voxel **x** in image **k** holds the probability that voxel **x** has label **k**. Likelihood images must be between 0 and 1. Any values outside of the range are interpreted as the voxel being excluded from the voting. These voxels will be assigned label 0 in the output.

The problem is encoded in the form of energy minimization, consisting of a data term and a regularization term. The data term encodes the cost (penalty) associated with assigning the voxel **x** the label **k**. The parameter **mode** describes how likelihood images are mapped to the cost. 

* `VOTES_AGAINST` or `VA`. This mode is useful when the command is being used to combine several multi-label segmentations into a single one. Each likelihood image is assumed to be the proportion of segmentations that assign label *k* to voxel *x*. The data term equals to the error associated to assining the voxel *k* label *x*. This error is calculated as the sum of the likelihoods for all labels at *x* minus the likelihood for *k* at *x*. Note that the likelihoods do not have to add up to one, which may be interpreted as missing data for some voxels. 

* `LOG_LIKELIHOOD` or `LL`. The cost for label *k* at voxel *x* is the logarithm of the k-th likelihood image at *x*. This will assign infinite cost when the likelihood is zero. 

The regularlization term is encoded as **lambda** times the total number of neighboring voxels inside the mask (non-excluded region of the image) that have different labels. 

The optmization problem is solved using the Alpha-Expansion graph cut algorithm. Users of this functionality should cite the following papers. 

1. Yuri Boykov, Olga Veksler, Ramin Zabih, *Efficient Approximate Energy Minimization via Graph Cuts*, IEEE transactions on PAMI, vol. 20, no. 12, p. 1222-1239, 2001. 

2. Vladimir Kolmogorov and Ramin Zabih, *What Energy Functions can be Minimized via Graph Cuts?*, IEEE transactions on PAMI, vol. 26, no. 2, p. 147-159, 2004.

3. Yuri Boykov and Vladimir Kolmogorov, *An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision*, IEEE transactions on PAMI, vol. 26, no. 9, p. 1124-1137, 2004. 

As noted in the open source implementation of the graph cuts algorithms distributed under the General Public License, "This software can be used only for research purposes, you should cite the aforementioned paper in any resulting publication.  If you wish to use this software (or the algorithms described in the aforementioned paper) for commercial purposes, you should be aware that there is a US patent: R. Zabih, Y. Boykov, O. Veksler, *System and method for fast approximate energy minimization via graph cuts*, United Stated Patent 6,744,923, June 1, 2004.

The example below illustrates applying **-vote-mrf** with a user-specified mask. Voxels outside of the mask are first remapped to NaN (not a number) and thus excluded from the MRF optimization and given 0 label.

    c3d lhood01.nii lhood02.nii lhood03.nii mask.nii -popas M \
        -foreach -push M -replace 0 NaN -times -endfor \
        -vote-mrf VA 0.2 -o result.nii

#### -voxreg, -voxelwise-regression: Regression between two images

Syntax: `-voxreg regression_order `

Perform regression between corresponding voxels in two images. This command takes two images as input, X and Y. This command finds parameters b\_0, b\_1, ..., b\_k, such that Y is best approximated by b\_0 + b\_1 * X + b\_2 * X^2 + ... + b_k * X^k. Here is an example of linear regression. 

    $ c3d Y.nii X.nii -voxreg 2 
    REGCOEFF[0] = 5.56935
    REGCOEFF[1] = 0.844024

    $ c3d Y.nii X.nii -scale 0.844024 -shift 5.56935 -voxreg 2
    REGCOEFF[0] = 0
    REGCOEFF[1] = 1

#### -wrap: Wrap (rotate) image 

Syntax: `-wrap <vector> `

Wrap image around one or more voxel dimensions. Wrapping is typically used to correct for MRI wrap-around artifacts. The vector must have integer components, possibly negative. For example, 

    c3d badmri.nii.gz -wrap 0x20x0 -o fixedmri.nii.gz

will wrap the image in the second voxel dimension by 20 voxels (i.e., voxel at 10x40x20 will me moved to the position 10x20x20). 


### Commands: Options and Parameters

Options change the behavior of commands that *appear later on the command line*. This is very important. Specifying options after the command will have no effect.

#### -background: Specify background intensity

Syntax: `-background <value> `

Sets the background intensity for interpolation and other operations where some default background value is needed. Default is 0.

#### -compress, -no-compress: Enable/disable compression for some image files

Syntax: `-compress` or `-no-compress`

Turns on compressing for image file formats that support it. For some file formats, like NIFTI (.nii), compression is enabled automatically when the filename includes the **.gz** extension. For other formats, like MetaImage, you need to specify **-compress** to enable compression. The following two commands save the image as compressed NIFTI and MetaImage files:

    c3d input.nii -o output.nii.gz
    c3d input.nii -compress -o output.mha

#### -interpolation: Set interpolation mode

Syntax: `-interpolation <NearestNeighbor|Linear|Cubic|Sinc|Gaussian> [param]`

Specifies the interpolation used with **-resample** and other commands. Default is **Linear**. Gaussian interpolation takes as the parameter the standard deviation of the Gaussian filter (e.g, 1mm). Gaussian interpolation is very similar in result to first smoothing an image with a Gaussian filter and then reslicing it with linear interpolation, but is more accurate and has less aliasing artifacts. It is also slower, and should only be used with small sigmas (a few voxels across). 

Shorthand 0 can be used for *NearestNeighbor*, 1 for *Linear* and 3 for *Cubic*. For example:

    c3d -int 3 test.nii -resample 200x200x200% -o cubic_supersample.nii


#### -noround, -round: Floating point rounding behavior

Syntax: `-noround` or `-round `

By default, **convert3d** will round floating point values when converting to an integer, short or byte image. This command specifies that rounding should not be used. Rounding is used to avoid numerical errors stemming from the internal floating point representation. 

    c3d image1.img -type short -noround image2.img


#### -pim, -percent-intensity-mode: Set behavior of % specifier

Syntax: `-pim Quantile | q | ForegroundQuantile | fq | Range | r`

This options changes how the percent sign (%) is interpreted when specifying intensity values. **Quantile (q)** means that `10%` describes the 10th percentile of all intensity values in the image (i.e., 10% of the voxels have lower intensity). **ForegroundQuantile (fq)** is similar, but voxels with background intensity (see **-background** option) are excluded from the percentile computation. **Range (r)** changes the meaning of percent sign from percentile to the range between the minimum and maximum of the image, and `0.1%` becomes equal to MIN + 0.1 (MAX - MIN). The default is **Quantile**. 

    $ c3d comp01.png -verbose -pim Quantile -verbose -threshold 75% inf 1 0 
    Quantile 0.75 maps to 18

    $ c3d comp01.png -verbose -pim ForegroundQuantile -verbose -threshold 75% inf 1 0 
    Foreground quantile 0.75 (over 37467 voxels) maps to 58

    $ c3d comp01.png -verbose -pim Range -verbose -threshold 75% inf 1 0 
    Intensity range spec 0.75 maps to 191.25

#### -rf-param-patch: Random Forest training patch size

Syntax: `-rf-param-patch <size_spec>`

Set the radius of the patch used to generate features for the RF classifier. By default this is zero, which means that just the intensity of each voxel is used as a feature. Setting this to non-zero values will result in neighboring intensities also being used as features, and can improve classification in presence of complex image texture. The patch size in each dimension is (2 * radius + 1). See **-rf-train** command for details.

    # Set patch size to 5x5x5 
    c3d ... -rf-param-patch 2x2x2 ... -rf-train myforest.rf

#### -rf-param-usexyz: Random Forest coordinate features

Syntax: `-rf-param-usexyz`

Use the coordinates of voxels as additional features. This allows some geometric relations between different labels to be learned. Equivalent to the corresponding ITK-SNAP option.

#### -rf-param-treedepth: Random Forest tree depth

Syntax: `-rf-param-treedepth <integer>`

Sets the depth of the trees in the classifier. Default value is 30. Deeper trees can learn on more complex data but require more time. 

#### -rf-param-ntrees: Random Forest forest size

Syntax: `-rf-param-ntrees <integer>`

Sets the number of trees in the forest. Default value is 50. Larger forests are more robust but more time to train and apply. 

#### -spm, -nospm: SPM compatibility in Analyze output

Syntax: `-spm` or `-nospm `

These options specify whether use the SPM extension to the Analyze (.hdr,.img) format. When this option is on, the origin field stored by SPM in the Analyze header will be correctly interpreted. When saving analyze files, the origin will be set correctly. The default is equivalent to the **-nospm** option. Best to avoid this issue altogether by using NIFTI and SPM5 or later.

    c3d -spm in.hdr out.img.gz


#### -type: Specify pixel type for image output

Syntax: `-type < char | uchar | short | ushort | int | uint | float | double > `

Specifies the pixel type for the output image. By default, images are written in floating point (**float**) format. The type does not affect how images are processed, only how they are saved. 

    c3d image1.img -type short image2.img

Some images require data in certain types. For example, to save PNG images, uchar or ushort type must be specified.

#### -verbose: Enable verbose output of commands

Syntax: `-verbose`

Commands entered after the **-verbose** command will print debugging information. This can be turned off with **-noverbose**.

### Parameter Specifications


#### Vectors

Certain commands, like **-smooth**, take vectors as parameters. Vectors can be specified in voxels and in spacing units (mm). The following are some examples of allowed specifications: 

    c3d ... -smooth 1.0x1.2x0.9mm ...   // Anisotropic, specified in mm (spacing units)
    c3d ... -smooth 1.0mm               // Isotropic, specified in mm (spacing units)
    c3d ... -smooth 1.0x1.2x0.9vox ...  // Anisotropic, specified in voxels
    c3d ... -smooth 3vox                // Isotropic, specified in voxels

#### Geometry

Geometry specifications are similar to the ones used in ImageMagick. They have the format **S1xS2xS3+O1+O2+O3**, where S1,S2,S3 give the size of the region and O1,O2,O3 give the offset. The offset can be omitted, and will default to one. 

#### Intensity Values

Intensity values can be specified as floating point values. Additionally, some commands support values **-inf** and **inf** for negative and positive infinity. Furthermore, intensity can be specified as a percentile of overall image intensity (computed in the last image on the **convert3d** stack). Here are some examples using the **-clip** command: 

    c3d ... -clip 100 400                         // Floating point lower and upper bounds
    c3d ... -clip -inf 400                        // Infinite lower bound
    c3d ... -clip -inf 90%                        // Upper bound specified as a percentile

With the **-percent-intensity-mode** (shorthand **-pim**) option, you have fine control over how the percent sign is interpreted when specifying intensities. By default, `10%` means 10th percentile of all intensities in the image. With **-pim ForegroundQuantile**, the quantile is computed only over voxels that are different from the background value (specified with **-background** option. This is useful if you need a percentile over a masked region of an image, etc. Lastly, with **-pim Range**, the percent sign changes meaning from percentile to the range between the minimum and maximum of the image, and `0.1%` becomes equal to MIN + 0.1 (MAX - MIN). 

### File Format Notes

**convert3d** supports all the file formats that ITK can read and write. In addition it supports some additional formats and format extensions. 

#### .hdr, .img, .img.gz 

Analyze files saved in SPM will by default be assigned a zero origin (as ITK does). But with the **-spm** option passed in before the filename, the origin will be correctly set. 

#### .cub 
VoxBo CUB files can be read and written 

#### .df3 
PovRay 3D density files can be written. They must be used with integral types (e.g., **-type ushort**) 

### Support for Multicomponent Images 

This section discusses images with multiple components, where each voxel contains more than one value. For example color (RGB) images and diffusion tensor images are multi-component images. C3D can read and write multi-component images, but at the present, it can only represent these images internally as separate single-component images. 

Use the **-mcs** flag to tell **convert3d** that when it encounters a multi-component image file, it should load each of the components as separate images. If you don't use this flag, only the first component of a multi-component image will be read in. 

    c3d -mcs ...

For example, suppose *rgb.mha* is a three-component color image. If you read it with the **-mcs** flag, **convert3d** will create three volumes, one for the red component, one for the green and one for the blue. You can use the **-foreach** construct to apply operations to the components. You can then save the components into a multi-component image with the **-omc** command. Here is an example: 

    c3d -mcs rgb.mha -foreach -scale 2 -shift 1 -endfor -omc 3 rgbscaled.mha

Here is another example, where we take an RGB image and extract red, green and blue channels. 

    c3d -verbose rgb.mha -o blue.nii.gz -pop -o green.nii.gz -pop -o red.nii.gz

The components are being saved in reverse order. That's because when the image is read, the red component is placed at the top of the stack, then the green component is placed on the top of the red component, and finally, the blue component is placed at the top of the green component. So the blue component is at the top of the stack. Since **convert3d** commands like **-o** apply to the image at the top of the stack, the blue image is saved first, then removed (popped off) the stack, followed by the green and blue components. 

With the advent of the **-oo** command (output multiple images), the command line becomes a little cleaner: 

  c3d -mcs rgb.mha -oo red.nii.gz green.nii.gz blue.nii.gz

or, if you want to name the component images as **comp00.nii.gz**, **comp01.nii.gz**, **comp02.nii.gz**: 

    c3d -mcs rgb.mha -oo comp%02d.nii.gz            

To create an RGB image from component images, we would use the following command: 

    c3d red.nii.gz green.nii.gz blue.nii.gz -omc 3 rgb.mha

### Commands with Variable Number of Parameters

Some more recently added commands in **c3d** take multiple parameters and apply them to multiple images on the stack. As a rule, the last parameter is applied to the topmost image on the stack, the second to last parameter is applied to second from the top image on the stack, and so on. 

For example, the \``\`-weighted-sum''' command computes a weighted sum of several images: 

    c3d A.nii B.nii C.nii D.nii -weighted-sum 0.1 0.2 0.3 0.4 -o W.nii

The result of this command will be 0.1 * A + 0.2 * B + 0.3 * C + 0.4 * D. The equivalent way to compute this expression is 

    c3d A.nii -scale 0.1 B.nii -scale 0.2 -add C.nii -scale 0.3 -add D.nii -scale 0.4 -add -o W.nii

Another example is the **-oo** command, which allows us to save several images at once. 

### How to Add a Command to C3D (For Software Developers)

This section describes the steps you must typically take to add a new command to C3D. Almost all commands are implemented in an adapter object. An adapter is simply a C++ class that inherits from `ConvertAdapter` and implements the `()` operator, where the work of the command is performed. There are many examples of adapters in the `adapters` directory. Some very simple commands, as well as options, stack manipulation, and flow control commands do not have adapters. We will discuss adding a command that is implemented by an adapter. 

1. To add a command to C3D, you should first create an adapter object. To simplify things, we provide a bash script to create a new adapter in the `adapters/generator` folder. Simply pass in the name of the generator and the parameters to the `operator ()` to the script. For example, if we are adding a command with a single parameter that is a vector, we would call the script as follows: 

        cd MYDIR/convert3d/adapters/generator
        bash runme.sh MyNewAdapter "RealVector x"
        cat MyNewAdapter.cxx                              <- Check that everything is OK
        cat MyNewAdapter.h                                 
        mv MyNewAdapter.* ..                              <- Copy to the adapter directory

2.  Next, you must write the code for the command. Typically, all you have to do now is to fill out the code for the `operator ()`. Please be sure to use the same convention as in other adapters. For example, every adapter should use the `*c->verbose` stream to report what it's doing and which image it's doing it to. Commands in C3D consume their arguments, so be sure to pop your images off the stack and push your result on the stack. Also be sure to do necessary error checking and fire off `ConvertException` if parameters or the state of the stack are invalid. 

3.  Next, you must add your new adapter to the `CMakeLists.txt` file. That's easy, just look for all the other adapters and insert yours in alphabetical order. 

4.  Next, edit `ConvertImageND.cxx`. This is the main application driver and command line processor. At the top of the file, add an `#include` statement for your new adapter, keeping things in alphabetical order. 

5.  Staying in `ConvertImageND.cxx`, insert your new command in the `usage()` function, again, keeping things in alphabetical order. 

6.  Still in `ConvertImageND.cxx`, find the section of the code where commands are implemented. Find where your command fits in alphabetical order, and insert a section of code corresponding to your command. Again, use existing commands as examples. Make sure to return a value - how many parameters your command has consumed from the command line. Otherwise, you will get a "fell through" error message. 

7.  Almost done. Edit the `utilities/bashcomp.sh` script and add a line for your command there. 

8.  Last step: edit the documentation file (c3d.md) to add your command's description

To summarize, here is a checklist for adding a new command 

1.  Generate adapter via bash script 
2.  Write code for `operator ()` 
3.  Add line in `CMakeLists.txt` 
4.  Add `#include` line in `ConvertImageND.cxx` 
5.  Add `usage()` line in `ConvertImageND.cxx` 
6.  Add code to parse command and call adapter in `ConvertImageND.cxx` 
7.  Add line to `bashcomp.sh` 
8.  Document the command on the wiki

 [1]: http://www.itksnap.org/pmwiki/pmwiki.php?n=Convert3D.Convert3D
 [2]: http://en.wikipedia.org/wiki/Reverse_Polish_notation
 [3]: http://c3d.cvs.sourceforge.net/viewvc/*checkout*/c3d/convert3d/utilities/bashcomp.sh
 [4]: http://www.insight-journal.org/browse/publication/640
 [5]: http://en.wikipedia.org/wiki/Morphological_image_processing
 [6]: http://www.fil.ion.ucl.ac.uk/spm/doc/books/hbf2/pdfs/Ch7.pdf
 [7]: http://www.itk.org/Doxygen/html/classitk_1_1BinaryFillholeImageFilter.html
 [8]: http://www.cppreference.com/wiki/c/io/printf
 [9]: http://www.itksnap.org/pmwiki/pmwiki.php?n=Convert3D.Point
 [10]: #ParameterSpecifications
 [11]: http://ieeexplore.ieee.org/xpls/abs_all.jsp?arnumber=1309714
 [12]: http://picsl.upenn.edu/ANTS/index.php
 [13]: http://www.imagemagick.org/script/command-line-options.php#colorspace
 [14]: http://www.itk.org/Doxygen/html/classitk_1_1CannyEdgeDetectionImageFilter.html
 [15]: http://www.itk.org/Doxygen/html/classitk_1_1RecursiveGaussianImageFilter.html
 [16]: http://www.insight-journal.org/browse/publication/721
