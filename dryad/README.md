# X-Ray CT scans of barley panicles and their individual seeds from the Composite Cross II experiment

## General information

### Authors

- **Erik J. Amézquita**, _Michigan State University_
- **Michelle Quigley**, _Michigan State University_
- **Tim Ophelders**, _TU Eindhoven_
- **Jacob B. Landis**, _Cornell University_
- **Daniel Koenig**, _University of California Riverside_
- **Elizabeth Munch**, _Michigan State University_
- **Daniel H. Chitwood**, _Michigan State University_

### To whom correspondence should be addressed:

**Erik J. Amézquita**
428 S Shaw Ln
Engineering Building Rm 1515
East Lansing, MI 48824
USA
amezqui3@msu.edu

**Daniel H. Chitwood**
1066 Bogue St
East Lansing, MI 48824
USA
+1 517 353 0462
chitwoo9@msu.edu

### Date and geographic location of data collection

In November of 2016, seeds from each accession were stratified at 4C on wet paper towels for a week, and germinated on the bench at room temperature at University of California, Riverside. Four day old seedlings were transferred into pots in triplicate and arranged in a completely randomized design in a greenhouse. Day length was extended throughout the experiment using artificial lighting ---minimum 16h light / 8h dark. After the plants reached maturity and dried, a single spike was collected from each replicate for scanning at Michigan State University throughout Fall 2018.

### Keywords

- Barley (Hordeum vulgare)
- inflorescence
- X-ray computed tomography (CT)
- Plant morphology
- 3D morphology

### Related content

- Amézquita EJ, Quigley M, Ophelders T, Landis JB, Koenig D, Munch E, Chitwood DH (2021) "Measuring hidden phenotype: Quantifying the shape of barley seeds using the Euler Characteristic Transform" [Preprint DOI](https://doi.org/10.1101/2021.03.27.437348). To appear in _inSilico Plants_.

- [Demeter scripts](https://github.com/amezqui3/demeter)

### License

CC0 1.0 Universal (CC0 1.0)
Public Domain Dedication 

[See summary here](https://creativecommons.org/publicdomain/zero/1.0/).

### Acknowledgements

DC is supported by the USDA National Institute of Food and Agriculture, and by Michigan State University AgBioResearch. The work of EM is supported in part by the National Science Foundation through grants CCF-1907591 and CCF-2106578. JL was supported by the NSF Plant Genome Postdoctoral Fellowship 1711807. DK is supported by an award from the National Science Foundation Plant Genome Research Program (IOS-2046256) and funding from the USDA NIFA (CA-R-BPS-5154-H).

=========

## Data and file overview

### Overview

We selected 28 barley accessions with diverse spike morphologies and geographical origins for our analysis ([Harlan and Martini 1929, 1936, 1940](https://doi.org/10.2134/agronj1929.00021962002100040014x)). In November of 2016, seeds from each accession were stratified at 4C on wet paper towels for a week, and germinated on the bench at room temperature. Four day old seedlings were transferred into pots in triplicate and arranged in a completely randomized design in a greenhouse. Day length was extended throughout the experiment using artificial lighting ---minimum 16h light / 8h dark. After the plants reached maturity and dried, a single spike was collected from each replicate for scanning at Michigan State University.
The scans were produced using the North Star Imaging X3000 system and the included efX software, with 720 projections per scan, with 3 frames averaged per projection. The data was obtained in continuous mode. The X‐ray source was set to a voltage of 75 kV, current of 100 &mu;A, and focal spot size of 7.5&mu;m. The 3D reconstruction of the spikes was computed with the efX-CT software, obtaining a final voxel size of 127 microns. The intensity values for all raw reconstructions was standardized as a first step to guarantee that the air and the barley material had the same density values across all scans. Next, the air and debris were thresholded out, and awns digitally.
Finally, the seed coat of the caryopses was digitally removed, leaving only the embryo and endosperm due to their high water content. We did not have enough resolution in the raw scans to distinguish clearly the endosperm from the embryo. Hereafter, we will refer to these embryo-endosperm unions simply as seeds. Due to the large volume of data, we used an in-house scipy-based python script to automate the image processing pipeline for all panicles and grains.

### File description

The whole dataset is split into four (4) clusters. Read the individual README files within each cluster, and see Amézquita et al (2021) for more information.

- `demeter`: collection of in-house python-based and R-based scripts to process the images, and later compute tradtional and topological shape descriptors. Tutorials on how to use the different functions are provided as commented jupyter notebooks. Created mainly spring 

- `dimension_reduction`: collection of CSV files with the topological shape descriptors for each of the 3121 founder seeds. These descriptors have been reduced in dimension using independently [KPCA](http://dx.doi.org/10.1162/089976698300017467) and [UMAP](http://arxiv.org/abs/1802.03426).

- `shape_descriptors`: collection of CSV files with the tradtional and topological shape descriptors for all the seeds. For ease of computation, part of the information is redundant. Shape descriptors are provided as CSVs for the 3121 founder seeds and for the 37 881 total seed collection.

- `xray_ct_scans`: 774 individual barley panicle and their corresponding 37 881 clean, individual seeds. All the scans are provided as single 3D 8-bit TIFF files. Three raw X-ray CT scans containing 4 barley panicles each are included as well.

- `corrected_metadata.csv`: CSV with scanning information of the 774 barley panicles as part of 224 raw scans.

- `LICENSE`: raw text file with Open Data Commons Attribution License (ODC-By) details

- `LICENSE_summary`: raw text with ODC-By human-readable summary.

- `README.md`: This file. Markdown format. Raw text.
