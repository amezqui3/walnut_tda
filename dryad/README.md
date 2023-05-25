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

The whole dataset is split into three (3) folders, plus additional metadata in CSV format. Read the individual README files within each folder for more details, and see Amézquita _et al._ (forthcoming) for more information.

- `merged`: Collection of summaries in CSV format. These summaries include the volume of all the main tissues of all citrus fruit, ellipsoid semi-axes lengths that model the overal fruit shape, fruit sphericity for all samples, among other measurments.

- `oil`: Collection of 8-bit TIFF files containing all the segmented oil glands for each citrus. Additional numbers are extracted related to the oil glands. Based on the oil gland centers, we compute the ellipsoid surface that approximates the best the overall fruit shape. Parameters of this best-fit ellipsoid are provided. See the overview section above for more details.

- `tissue`: Separate 8-bit TIFF files for separate tissues. All images have the same shape as the original whole fruit image. These are central column, endocarp, rind, and exocarp for each fruit. A rotation matrix to align the fruit upright is provided.

- `CRC_citrus_scan_technique.csv`: CSV with resolution of the scans.

- `CRC_citrus_scanned.csv`: Metadata provided by the Walnut Improvement Program at the University of California Davis. 

- `LICENSE`: raw text file with CC0 License details

- `LICENSE_summary`: raw text with CC0 human-readable summary.

- `README.md`: This file. Markdown format. Raw text.
