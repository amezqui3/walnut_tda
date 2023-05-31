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

Collection of 1-row, headless CSVs with several phenotypes of the walnut and its shell, kernel, and packing tissues.

See the jupyter notebook(s) for more details:

```
- https://github.com/amezqui3/walnut_tda/blob/main/jupyter/05_shell_thickness.ipynb
- https://github.com/amezqui3/walnut_tda/blob/main/jupyter/07_traditional_phenotyping.ipynb
- https://github.com/amezqui3/walnut_tda/blob/main/jupyter/08_phenotype_wrangling.ipynb
```

Subdirectories are named as `<Locate>_<Row>_<Tree>` to idenfify individual walnuts from different trees. See the main README file for more details.

- `_trad.csv`: Single-line, headless CSV with 48 entries.
    - UCACCSD identifier. See README file for more details
    - `<Locate>_<Row>_<Tree>` Identifier. See README file for more details.
    - Individual label
    - Nut Length along the proximal-distal axis [mm]
    - Nut Height along the medial-lateral axis [mm]
    - Nut Width along the adaxial-abaxial axis [mm]
    - Nut Volume [mm^3]
    - Nut VA3D sphericity index [1 - infinity, where 1 indicates perfect sphericity]
    - Nut Feret Ratio [same as above]
    - Nut Surface Area [mm^2]
    - Nut Sphericity, Wadell index [0 - 1, where 1 indicates perfect sphericity and 0 not spherical at all]
    - Surface area of the convex hull of the nut [mm^2]
    - Volume of the convex hull [mm^3]
    - Ratio between the surface area of the convex hull and the surface area of the nut [0 - 1, where 1 indicates perfect convexity and 0 not convex at all]
    - Ratio between the volume of the nut and the volume of its convex hull [same as above]
    - Krumbein Index of sphericity [same as above]
    - Sneed Index of sphericity [same as above]
    - Air Volume [mm^3]
    - Kernel Volume [mm^3]
    - Shell Volume [mm^3]
    - Packing tissue Volume [mm^3]
    - Air Volume Ratio, how much of the total volume is air [%]
    - Kernel Volume Ratio, how much of the total volume is kernel [%]
    - Shell Volume Ratio [%]
    - Packing tissue Volume Ratio [%]
    - Shell Rugosity index [1 - infinity, where 1 indicates perfect smoothness] 
    - Average thickness of the shell across the whole walnut [mm]
    - Protruding Shell Volume Ratio, how much of the total volume of the shell protrudes into the interior of the walnut [%]
    - Protruding Shell Volume [mm^3]
    - Kernel Length along the proximal-distal axis [mm]
    - Kernel Height along the medial-lateral axis [mm]
    - Kernel Width along the adaxial-abaxial axis [mm]
    - Kernel Surface Area [mm^2]
    - Volume of the convex hull of the kernel [mm^3]
    - Surface area of the convex hull of the kernel [mm^2]
    - Ratio between the surface area of the convex hull and the surface area of the kernel [0 - 1, where 1 indicates perfect convexity and 0 not convex at all]
    - Ratio between the volume of the kernel and the volume of its convex hull [0 - 1, where 1 indicates perfect convexity and 0 not convex at all]
    - Surface area of the main cavity surrounded by the kernel at its distal end [mm^2]
    - Volume of the main cavity surrounded by the kernel at its distal end [mm^3]
    - Length along the proximal-distal axis of the main arch that connects the two hemispheres of a walnut [mm]
    - Length along the proximal-distal axis of the main cavity surrounded by the kernel at its distal end [mm]
    - Wadell's index of sphericity for the kernel [0 - 1, where 1 indicates perfect sphericity and 0 not convex at all]
    - Ratio between the surface area of the main cavity and the surface area of the kernel [0 - 1]
    - Ratio between the volume of the main cavity and the volume of the convex hull of the kernel [0 - 1]
    - Ratio between the average raw density of the kernel and the average raw density of the shell
    - Ratio between the average raw density of the packing tissue and the average raw density of the shell
    - Ratio between the average raw density of the packing tissue and the average raw density of the kernel
    - Feret ratio of the kernel [1 - infinity, where 1 indicates perfect sphericity]
    
