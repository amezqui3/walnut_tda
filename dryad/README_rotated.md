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

All plant materials represent walnut breeding lines, germplasm, and cultivars maintained by the Walnut Improvement Program at the University of California, Davis. A total of 150 walnuts accessions were harvested into mesh bags at hull split, oven-dried overnight at 95F, and then air-dried for several weeks before moving into cold storage at 35F. 5 to 16 individuals were selected for each accession, for a total of 1301 individual walnuts to be scanned at Michigan State University (Table S1). The walnuts were scanned in 171 batches. The scans were produced using the the North Star X3000 system and the included efX-DR software. The X-ray source was set at 75~kV and 100~\micro{A}, with 720 projections per scan, at 3 frames per second and with 3 frames averaged per projection. The data was obtained in continuous mode. The 3D X-ray CT reconstruction was computed with the efX-CT software, obtaining voxel-based images with voxel size of 75.9 microns.

### File description

Rotation and centering matrices and related data for the walnuts. Many phenotypes depend on having all the walnuts aligned the same way. In this case, we make sure that the X axis corresponds to the proximal-distal axis, with the tip of the walnut on the postive side. The seal of the walnut lies on the XZ plane. The seal is parallel to the Y axis while perpendicular to the Z one. 

See the jupyter notebook(s) for more details:

```
- https://github.com/amezqui3/walnut_tda/blob/main/jupyter/03_nut_alignment.ipynb
```

Subdirectories are named as `<Locate>_<Row>_<Tree>` to idenfify individual walnuts from different trees. See the main README file for more details.

- `_rotation.csv`: 21x3 array with the following information by row. All measurements in milimeters. Using python notation, if this CSV is loaded as a numpy array `rdata`, then
    - `rdata[0]`: X,Y,Z coordinates of the center of mass of the walnut
    - `rdata[1:4]`: Rotation matrix on the X axis so that the seal of the walnut is aligned with the Y axis and perpendicular to the Z one.
    - `rdata[4:7]`: Rotation matrix on the Y axis so that the tip of the walnut is on the X axis
    - `rdata[7:10]`: Rotation matrix on the Z axis so that the tip of the walnut is on the X axis
    - `rdata[10:13]`: Final rotation matrix. If used, the 3 prior matrices should be ignored.
    - `rdata[13]`: Bools on wheter the X,Y,Z was flipped. Only `rdata[13,0]` is of importance, as the only possible flip that could have happened is having the tip of the walnut on the positive side of the X-axis rather than on the negative one.
    - `rdata[14]`: XYZ Coordinates of the tip of the walnut on the original scan. Identifying such tip guides much of the later rotations.
    - `rdata[15]`: XYZ coordinates of the tip of the walnut after this has been rotated.
    - `rdata[16]`: Most negative X coordinate, most negative Y coordinate, and most negative Z coordinate.
    - `rdata[17]`: Most positive X, Y, Z coordinates
    - `rdata[18]`: Feret (caliper) diameter along the X, Y, Z axes
    - `rdata[19]`: Surface area of the convex hull, volume of the convex hull, 0
    
- `_original.png`: 2D projection of the original, clean walnut

- `+rotXYZ.png`: 2D projection of the same walnut after rotation.
