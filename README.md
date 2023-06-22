# ETH3D Format Loader Example

Demonstrates how to load camera extrinsics and intrinsics for datasets of the [ETH3D benchmark](https://www.eth3d.net/).
The file format is the [human-readable COLMAP text format](http://colmap.github.io/format.html#text-format).

Usage:

```
ETH3DFormatLoader <path_to_cameras_txt> <path_to_images_txt>
```

If you use the dataset for research, please cite our paper:

T. Schöps, J. L. Schönberger, S. Galliani, T. Sattler, K. Schindler, M. Pollefeys, A. Geiger, "A Multi-View Stereo Benchmark with High-Resolution Images and Multi-Camera Videos", Conference on Computer Vision and Pattern Recognition (CVPR), 2017. \[[Bibtex](https://www.eth3d.net/data/schoeps2017cvpr.bib)\]\[[PDF](https://www.eth3d.net/data/schoeps2017cvpr.pdf)\]\[[Supplementary](https://www.eth3d.net/data/schoeps2017cvpr-supp.pdf)\]