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

We selected 28 barley accessions with diverse spike morphologies and geographical origins for our analysis ([Harlan and Martini 1929, 1936, 1940](https://doi.org/10.2134/agronj1929.00021962002100040014x)). 

In November of 2016, seeds from each accession were stratified at 4C on wet paper towels for a week, and germinated on the bench at room temperature. Four day old seedlings were transferred into pots in triplicate and arranged in a completely randomized design in a greenhouse. Day length was extended throughout the experiment using artificial lighting ---minimum 16h light / 8h dark. After the plants reached maturity and dried, a single spike was collected from each replicate for scanning at Michigan State University. The scans were produced using the North Star Imaging X3000 system and the included efX software, with 720 projections per scan, with 3 frames averaged per projection. The data was obtained in continuous mode. The X‐ray source was set to a voltage of 75 kV, current of 100 &mu;A, and focal spot size of 7.5&mu;m. The 3D reconstruction of the spikes was computed with the efX-CT software, obtaining a final voxel size of 127 microns. The intensity values for all raw reconstructions was standardized as a first step to guarantee that the air and the barley material had the same density values across all scans. Next, the air and debris were thresholded out, and awns digitally. 

Finally, the seed coat of the caryopses was digitally removed, leaving only the embryo and endosperm due to their high water content. We did not have enough resolution in the raw scans to distinguish clearly the endosperm from the embryo. Hereafter, we will refer to these embryo-endosperm unions simply as seeds. Due to the large volume of data, we used an in-house scipy-based python script to automate the image processing pipeline for all panicles and grains.

### File description

This is a collection of CSVs with traditional and topological shape descriptors for each file. The topological shape descriptors are based off the [Euler Characteristic Transform (ECT)](http://dx.doi.org/10.1093/imaiai/iau011). See Amézquita et al (2021) for more details. 

To ease the analysis pipeline, information in these files contains several redundancies across files.

- `founder_seeds`: Shape descriptors for all 3121 clean, founder seeds.
   - `combined_d<d>_T<T>.csv`: Metadata, traditional, and topological shape descriptors (in that order) for every founder seed.
   - The metadata (9 columns) is imported from `corrected_metadata.csv`
   - The traditional shape descriptors (11 columns)
   - The values of the ECT computed for `d` directions and `T` thresholds, resulting in `d x T` columns.
   - `T = 64` in this case.
   - Observe that the ECT values for `T = 32` can be easily obtained by just considering every other column of the topological descriptors.
   - Similarly, `T = 16` corresponds to take one topological shape descriptor column out of every 4.
   - And so on.

- `all_seeds`: Shape descriptors for all 37 881 files.
    - `all_seeds_traditional_raw.csv`: Traditional shape descriptors for the initial collection of seeds, before manual examination and cleaning.
    - `all_seeds_traditional_clean.csv`: Traditional shape descriptors for the final, clean collection of 37 881 seeds.
    - `all_seeds_combined_d158_T64`: Same idea as in the `founder_seeds` case, except now we consider the whole seed collection.
