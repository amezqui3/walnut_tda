# The shape of kernels and cracks, in a nutshell

## General information

### Authors

- **Erik J. Amézquita**, _Michigan State University_
- **Michelle Quigley**, _Michigan State University_
- **Patrick J. Brown**, _University of California Davis_
- **Elizabeth Munch**, _Michigan State University_
- **Daniel H. Chitwood**, _Michigan State University_

### To whom correspondence should be addressed:

**Daniel H. Chitwood**
1066 Bogue St
East Lansing, MI 48824
USA
+1 517 353 0462
chitwoo9@msu.edu

### Date and geographic location of data collection

Data collected in December 2019 at Michigan State University. Walnut samples provided by the [Walnut Improvement Program](https://fruitsandnuts.ucdavis.edu/collaborators/california-walnut-board/reports) at University of California Davis

### Keywords

- Persian Walnut (Juglans regia)
- Walnut
- X-ray computed tomography (CT)
- Plant morphology
- 3D morphology

### Related content

- Amézquita EJ, Quigley M, Brown PJ, Munch E, Chitwood DH "The shape of kernels and cracks, in a nutshell" In preparation at the time of writing this README.

- [Walnut scripts](https://github.com/amezqui3/walnut_tda)

### License

CC0 1.0 Universal (CC0 1.0)
Public Domain Dedication 

[See summary here](https://creativecommons.org/publicdomain/zero/1.0/).

### Acknowledgements

Daniel Chitwood is supported by the USDA National Institute of Food and Agriculture, and by Michigan State University AgBioResearch. The work of Elizabeth Munch is supported in part by the National Science Foundation through grants CCF-1907591, CCF-2106578, and CCF-2142713.

=========

## Data and file overview

### Overview

We explore the shape of walnut shells and kernels. Here, we study the shape of walnut fruits based on the X-ray CT 3D reconstruction of 1264 different samples comprising 150 accessions maintained by the Walnut Improvement Program at the University of California Davis. We exploit the nondestructiveness of X-rays to digitally segment and measure the 4 main tissues of interest for each walnut, namely shell, kernel, packing tissue, and sealed air. From these we extract a total of 38 size- and shape-specific descriptors, many of them unexplored in the current literature. We focus on several allometric relationships of interest, from which we draw theoretical upper and lower bounds of possible walnut and kernel sizes. We then study correlations and variations of these morphological descriptors with qualitative data corresponding to traits of commercial interest like ease of kernel removal and shell strength.

All plant materials represent walnut breeding lines, germplasm, and cultivars maintained by the Walnut Improvement Program at the University of California, Davis. A total of 150 walnuts accessions were harvested into mesh bags at hull split, oven-dried overnight at 95F, and then air-dried for several weeks before moving into cold storage at 35F. 5 to 16 individuals were selected for each accession, for a total of 1301 individual walnuts to be scanned at Michigan State University (Table S1). The walnuts were scanned in 171 batches. The scans were produced using the the North Star X3000 system and the included efX-DR software. The X-ray source was set at 75~kV and 100~\micro{A}, with 720 projections per scan, at 3 frames per second and with 3 frames averaged per projection. The data was obtained in continuous mode. The 3D X-ray CT reconstruction was computed with the efX-CT software, obtaining voxel-based images with voxel size of 75.9~{\micro}m.

### File description

Collection of 8-bit TIFF files containing all the segmented tissues of individual walnuts as voxel-based X-ray CT scans. These tissue are the kernel, shell, and packing tissues. Air contained in the walnut is also segmented. For all the tissue identified as shell, we also idenfied which voxels protruded inside the walnut cavity.

Individual tissues where roughly separated at first by a combination of thresholding and morphological operators (erosions, dilations). These operations where based on location, thickness, and density of voxels. These produced rough estimates of where each tissue lied. The paramters used were the same for all the walnuts. Second, this estimates were completed with watershed segmentation.

A visual inspection was later peformed to assess the segmentation correctness and manual recalibration of image processing parameters was done for a handful of walnuts.

These scans are already clean and standardized. They have no external air (only the voxels corresponding to air within the walnut cavity report nonzero values). All the planes containing exclusively zeros were removed in order to reduce image dimensions. The density values for shell, kernel, and interior air are roughly the same across samples. This standarization helps later when separating tissues with a watershed approach. 

See the jupyter notebook(s) for more details:

```
- https://github.com/amezqui3/walnut_tda/blob/main/jupyter/02_watershed_segmentation.ipynb
- https://github.com/amezqui3/walnut_tda/blob/main/jupyter/04_interior_shell.ipynb
```

Subdirectories are named as `<Locate>_<Row>_<Tree>` to idenfify individual walnuts from different trees. See the main README file for more details.

Inside each subdirectory, TIFF files are located, indenfied with suffixes `001`, `002`, etc. whenever they correspond to the same individual walnut.

- `_air.tif`: The air contained inside the walnut
- `_meat.tif`: The kernel contained in the walnut
- `_shell.tif`: The shell of the walnut, including highly dense, connected voxels that protrude inside the walnut cavity.
- `_vein.tif`: The packing tissue inside the walnut.
- `_protrusion.tif`: 4-color-coded version of the `_shell.tif` scan.
    - `0`: Non-shell material
    - `1`: Shell voxel close to the border of the walnut. Corresponds the intuitive idea of what the shell should be.
    - `2`: Shell-like voxel that is part of a strange bulge that goes a bit into the walnut. We are still not 100% sure what this is. The voxel has the same density as the shell close to the border. And it is connected to such border by similarly-dense voxels.
    - `3`: Shell voxel that clearly protrudes inside the walnut from the proximal end.

Inside each subdirectory, extra folders are located with additional information for each individual walnut scan.

- `diagnostic/`: Quick visual assessments of 2D slices of walnut scan across the X,Y,Z axes. Columns show: 
    - initial, clean, invidiual scan of the whole walnut
    - Color-coded tissues after the watershed segmentation
    - Only the kernel
    - Only the shell
    - Only the packing tissue
    
- `protrusion/`: Quick visual assessments of 2D projections and slices of the shell and the identification of protruding shell material.
    - `_cavity.jpg`: 2D projections of the `_protrusion.tif` scan. External shell is in purple. Shell-like material that bulges into the walnut is in orange. Shell-like material that protrudes from the proximal end is yellow.
    - `_coords.jpg`: Locates the tips of all the protruding bits of shell-like tissue. Bits that protrude from the proximal end are in blue. The rest are red.
    - `_detection.jpg`: All the shell and shell-like material 2D slices. Colored in pale yellow. Areas of where protrusion from the proximal end might be going on are in purple.
    - `_detection.jpg`: All the shell and shell-like material 2D slices. External shell in purple. Shell protruding from the proximal end in yellow. Shell bulging into the walnut from elsewhere in orange.
    
- `snaps/`: Quick visual assessment of 2D projections of the kernel, air, shell, and packing tissues. See original README for more details. Name of the JPGs corresponds with the TIFF scans. 
