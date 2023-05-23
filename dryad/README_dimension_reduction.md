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

This is a collection of CSVs of topological shape descriptors reduced in dimension for each of the 3121 founder seeds. The topological shape descriptors are based off the [Euler Characteristic Transform (ECT)](http://dx.doi.org/10.1093/imaiai/iau011). These descriptors have been reduced in dimension using independently [KPCA](http://dx.doi.org/10.1162/089976698300017467) and [UMAP](http://arxiv.org/abs/1802.03426). See Amézquita et al (2021) for more details. 

To ease the analysis pipeline, information in these files contains several redundancies across files.

- `kpca`: Traditional and Kernel PCA dimension-reduced topological shape descriptors for each founder seed.
    - The first 11 columns correspond to traditional shape descriptors
    - The rest of columns correspond to dimension-reduced topological shape descriptors
    - All descriptors have been centered and scaled to variance 1.
    - KPCA was computed with Laplacian kernel and `sigma = 1`. Refer to the corresponding Jupyter notebook in the `demeter` respository for more details.
    - `kpca_normalized_size_combined_<d>_<T>_laplacedot_<dim>_founders.csv`
    - `d` is the number of directions used to compute the ECT
    - `T` is the number of thresholds used to compute the ECT
    - `dim` is the number of dimensions that the ECT was reduced to
    - Just like traditional PCA, reducing the ECT to just 2 dimensions is equivalent to first reduce the ECT to `dim = 24` dimensions, and then just consider its first 2 principal components.
    
    - `kpca_normalized_size_traditional_158_8_vanilladot_0_founders.csv`: Linear PCA applied to just the traditional shape descriptors.
    
- `umap`: Topological shape descriptors reduced in dimension with unsupervised UMAP. All founder seeds.
    - `umap_<d>_<T>_<n>_<min>_<dim>_<metric>_unsupervised.csv`
    - `d` is the number of directions used to compute the ECT
    - `T` is the number of thresholds used to compute the ECT
    - `n` is the number of neighbors when computing UMAP
    - `min` is the minimum distance when computing UMAP
    - `metric` is the metric used to compute UMAP
    - `dim` is the number of dimensions that the ECT was reduced to
    
    
