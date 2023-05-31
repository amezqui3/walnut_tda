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

All plant materials represent walnut breeding lines, germplasm, and cultivars maintained by the Walnut Improvement Program at the University of California, Davis. A total of 150 walnuts accessions were harvested into mesh bags at hull split, oven-dried overnight at 95F, and then air-dried for several weeks before moving into cold storage at 35F. 5 to 16 individuals were selected for each accession, for a total of 1301 individual walnuts to be scanned at Michigan State University (Table S1). The walnuts were scanned in 171 batches. The scans were produced using the the North Star X3000 system and the included efX-DR software. The X-ray source was set at 75~kV and 100~\micro{A}, with 720 projections per scan, at 3 frames per second and with 3 frames averaged per projection. The data was obtained in continuous mode. The 3D X-ray CT reconstruction was computed with the efX-CT software, obtaining voxel-based images with voxel size of 75.9 microns.

### File description

Collection of CSVs with several phenotypes of the walnut kernel. We attempted to compute a rough approximation of the main cavity at the distal end of the kernel. We tried to gauge its size, both in terms of surface area and volume, both in absolute and relative terms. We emphasize that this is just a rough attempt and a more robust approach should be planned in this direction.

See the jupyter notebook(s) for more details:

```
- https://github.com/amezqui3/walnut_tda/blob/main/jupyter/06_kernel_traditional.ipynb
```

Subdirectories are named as `<Locate>_<Row>_<Tree>` to idenfify individual walnuts from different trees. See the main README file for more details.

- `_kernel.csv`: Single-line, headless CSV with 15 entries.
    - `<Locate>_<Row>_<Tree>` Identifier. See README file for more details.
    - Individual label
    - Surface area of the kernel [mm^2]
    - Volume of the kernel [mm^3]
    - Surface area of the kernel's convex hull [mm^2]
    - Volume of the convex hull [mm^3]
    - Krumbein sphericity index [0 - 1, where 1 indicates perfect sphericity and 0 not spherical at all]
    - Corey sphericity index [same as above]
    - Sneed sphericity index [same as above]
    - Janke sphericity index [same as above]
    - Equancy sphericity index [same as above]
    - Surface area of the main cavity at the distal end [mm^2]
    - Volume of the main cavity at the distal end [mm^3]
    - Height of the main arch that connects both hemispheres of the walnut [mm]
    - Height of the main cavity at the distal end [mm]
    
- `_cavity.jpg`: Visual assessment of the main cavity at the distal end of the kernel. Several 2D slices are take along the proximal-distal axis. The estimated cavity is drawn in red. The largest circle that can fit inside it is drawn in black.

- `_kernel.jpg`: Visual assessment of kernel, after it has been oriented. The length of the main arch between hemispheres is highlighted in red. The height of the main cavity at the distal end is marked in blue.
