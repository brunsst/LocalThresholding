# LocalThresholding

![alt tag](https://github.com/brunsst/LocalThresholding/TOC-Figure.png)

By default traverses a 3D image stack and calculates sparse Otsu thresholds every 50 voxels for a spherical region of 200 voxels radius.
Thresholds between the support points are interpolated linearly.
Default settings expect a 

### Compilation

Required libraries are LibTiff and OpenMP. The source code should compile on most Linux distributions by running: 

***<p align="center"> sh build_localthresholding.sh </p>***

which will provide an executable *locthresh* in the same directory.
