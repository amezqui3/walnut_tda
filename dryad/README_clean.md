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

Collection of 8-bit TIFF files containing all the individual walnuts as voxel-based X-ray CT scans. These scans are already clean and standardized. They have no external air (only the voxels corresponding to air within the walnut cavity report nonzero values). All the planes containing exclusively zeros were removed in order to reduce image dimensions.

The density values for shell, kernel, and interior air are roughly the same across samples. This standarization helps later when separating tissues with a watershed approach. 

See the jupyter notebook(s) for more details:

```
- https://github.com/amezqui3/walnut_tda/blob/main/jupyter/01_density_normalization.ipynb
```

Subdirectories are named as `<Locate>_<Row>_<Tree>` to idenfify individual walnuts from different trees. See the main README file for more details.

Inside each subdirectory, individual TIFF files are located, indenfied with suffixes `001`, `002`, etc.

Inside each subdirectory, a `normalization` folder is located with additional information for each individual walnut scan.

- `clean_zeros.csv`: 6-dimensional vector. This helps slicing the original raw image when removing planes containing exclusively zeros.

     `img = img[cero[1]:cero[4], cero[2]:cero[5], cero[0]:cero[3]]`
     
- `normalization_coefficients.csv`: 2-dimensional vector `[b,m]`. It is used to recalibrate the intensity values of the voxels of the raw scan 

     `new_density = m * raw_density + b`
     
- `walnut_shell.jpg`: 2D projection of the max intensity values of the whole walnut scan. Usually these max intensity voxel corresponds to one in the shell. See main README for more details.

- `whole_walnut.jpg`: 2D projection of the sum of intensity values of the whole walnut scan. See main README for more details.
