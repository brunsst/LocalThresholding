# LocalThresholding

![alt tag](https://github.com/brunsst/LocalThresholding/blob/main/TOC-Figure.png)

By default traverses a 3D image stack and calculates sparse Otsu thresholds every 50 voxels for a spherical region of 200 voxels radius.
Thresholds between the support points are interpolated linearly.
Default settings expect a tomographic reconstruction, i.e., zero values outside a cylindrical field of view are ignored.

### Compilation

Required libraries are LibTiff and OpenMP. The source code should compile on most Linux distributions by running: 

***<p align="center"> sh build_localthresholding.sh </p>***

which will provide an executable *locthresh* in the same directory.

### Demo

Running *locthresh* from the project directory with the *--demo* flag creates a folder */localseg/* with a segmentation of some sandgrains.

### Basic Arguments

Arguments can be passed via the command line. The following arguments are supported: 

| argument | value | explanation |
|--------|------------------|-----------|
| **-i** |/directory/with/tif/sequence/| (*mandatory*) input directory with tiff image sequence|
| **-o** |/directory/with/segmentation/| (*optional*) optional output directory. By default files are stored in */localseg/* at the level of input directory.|
| **-step** |integer| (*optional*) modifies the spacing of support points|
| **-radius** |integer| (*optional*) modifies the radius of the local (spherical) window|
| **-n_cpu** |integer| (*optional*) sets a maximum to the number of allowed parallel threads (default is 128).|
|**--subregion**|| (*optional*) image is not a cylindrical reconstruction. This includes zero values at the image boundaries.|
|**--demo**|| (*optional*) perform a segmentation on the files in the */demo/* directory.|
