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

- IPGRI (1994) _Descriptors for Walnut (_Juglans _spp.)_ International Plant Genetic Resources Institute, Rome, Italy. [URL](https://cgspace.cgiar.org/handle/10568/73159)

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

The whole dataset is split into five (5) folders, plus additional metadata in CSV format. Read the individual README files within each folder for more details, and see Amézquita _et al._ (forthcoming) for more information.

- `clean/`: Collection of 8-bit TIFF files containing all the individual walnuts as voxel-based X-ray CT scans. 
    - These scans are already clean and standardized. 
    - They have no external air (only the voxels corresponding to air within the walnut cavity report nonzero values). 
    - All the planes containing exclusively zeros were removed in order to reduce image dimensions.
    
- `kernel/`: Collection of CSVs with several phenotypes of the walnut kernel. 
    - We attempted to compute a rough approximation of the main cavity at the distal end of the kernel. 
    - We tried to gauge its size, both in terms of surface area and volume, both in absolute and relative terms. 
    - We emphasize that this is just a rough attempt and a more robust approach should be planned in this direction.
    
- `rotated/`: Collection of CSVs with rotation matrices, centering coordinates, and related data for the walnuts. 
    - Many phenotypes depend on having all the walnuts aligned the same way. 
    - In this case, we make sure that the X axis corresponds to the proximal-distal axis, with the tip of the walnut on the postive side. 
    - The seal of the walnut lies on the XZ plane. The seal is parallel to the Y axis while perpendicular to the Z one. 
    
- `traditional/`: Collection of 1-row, headless CSVs with several phenotypes of the walnut and its shell, kernel, and packing tissues.

- `all_phenos.csv`: A 721x20 table. Traits measured for breeding purposes. Most traits are measured on an ordinal scale following the criteria established by [(IPGRI, 1994)](https://cgspace.cgiar.org/handle/10568/73159). All these traits were measured on a per-tree basis by taking 10 different nuts. The columns are:
    - `UCACCSD`: identifier of the walnut cultivar according to the Walnut Improvement Program at the University of California Davis.
    - `YR`: year when walnuts were collected for evaluation of tree fruit traits
    - `LOCATE`: Location of the tree
    - `ROW`: Row where the tree is located
    - `TREE`: Exact tree within the row is located
    - `Shell Integrity`: Completeness of the shell
    - `Shell Texture`: Smoothness/roughness of the shell
    - `Shell Color`: Light or dark
    - `SEAL`: Open/very weak or very strong shell seal
    - `Shell Strength`: Paper or strong shell
    - `Shell thickness`: As measured near center of half shell [mm]
    - `Packing Tissue`: Thin/sparse or very thick packing tissue
    - `Kernel fill`: Poor or well filled
    - `Tip Shrivel`: Percentage of kernels exhibiting tip shrivel
    - `Minor Shrivel`: Percentage of kernels exhibiting minor shrivel
    - `Major Shrivel`: Percentage of kernels exhibiting major shrivel
    - `Plumpness`: Thin or plump kernel
    - `Ease of removal`: Easiness/difficulty to extract kernel halves
    - `Percent Kernel`: Percentage of the walnut's weight that corresponds to the kernel.
    - `Blank`: Percentage of nuts with no kernel

- `README.md`: This file. Markdown format. Raw text.

- `README_clean.md`: Copy of the README corresponding to the `clean` folder

- `README_kernel.md`: Copy of the README corresponding to the `kernel` folder

- `README_rotated.md`: Copy of the README corresponding to the `rotated` folder

- `README_traditional.md`: Copy of the README corresponding to the `traditional` folder

- `README_watershed.md`: Copy of the README corresponding to the `watershed` folder

