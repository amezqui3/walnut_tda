# X-Ray CT scans of barley panicles and their individual seeds from the Composite Cross II experiment

## General information

774 individual barley panicle and their corresponding 37 881 clean, individual seeds. All the scans are provided as single 3D 8-bit TIFF files. Three raw X-ray CT scans containing 4 barley panicles each are included as well. 

Scroll at the end for more content-specific details.

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

The barley spikes were scanned in batches of 4 most of the time. The panicles were proped up on a color-coded foam tray. 224 batches were scanned in total, producing 774 spikes. To keep track of the orientation and colors, the tray had a wooden letter attached to the red corner. Please refer to the included JPG images.

The scans were produced using the North Star Imaging X3000 system and the included efX software, with 720 projections per scan, with 3 frames averaged per projection. The data was obtained in continuous mode. The X‐ray source was set to a voltage of 75 kV, current of 100 &mu;A, and focal spot size of 7.5&mu;m. The 3D reconstruction of the spikes was computed with the efX-CT software, obtaining a final voxel size of 127 microns. The intensity values for all raw reconstructions was standardized as a first step to guarantee that the air and the barley material had the same density values across all scans. Next, the air and debris were thresholded out, and awns digitally.
Finally, the seed coat of the caryopses was digitally removed, leaving only the embryo and endosperm due to their high water content. We did not have enough resolution in the raw scans to distinguish clearly the endosperm from the embryo. Hereafter, we will refer to these embryo-endosperm unions simply as seeds. Due to the large volume of data, we used an in-house scipy-based python script to automate the image processing pipeline for all panicles and grains.

### File description

- `barley_8bit_tiff.jpg`: Image of the efX sofware when selecting 8-bit TIFF resolution.
- `barley_lab_composition.jpg`: Panel of the scanning pipeline. The X3000 NSI system, the foam tray to prop up the barley spikes and scan them by 4 at a time, screen shot of 3D reconstruction without air and foam.
- `barley_xray_setup.jpg`: Panicles before scanning.
- `foam_tray_*.jpg`: Closeup of the foam trays used.

- `corrected_metadata.csv`: CSV with scanning information of the 774 barley panicles as part of 224 raw scans.

- `LICENSE`: raw text file with Open Data Commons Attribution License (ODC-By) details
- `LICENSE_summary`: raw text with ODC-By human-readable summary.
- `README.md`: This file. Markdown format. Raw text.

- `raw_sample`: 3 raw TIFF scans. The scan contains no nonzero voxels. The air, foam, and other debris must be removed before further analysis.

- `raw_sample_auxiliary`: Auxiliary files to normalize the densities accross the 3 scans. See the corresponding density normalization Jupyter notebook in the `demeter` code respository.

- `spikes`: 3D 8-bit TIFF files. Contains 774 spikes and 224 wood markers. 
    - `S001` ... `S224`: spikes corresponding to each of the 224 batches.
    - Most of the batches contain 4 panicles. Refer to `corrected_metadata.csv` for more details.
    - Filenames of each scan: `S<scan>_l<label>_<color>_x<xcoord>_y<ycoord>_z<zcoord>.tif`
    - Scan number goes from `001` to `224`
    - Label number goes from `0` to `N`, where `N` is the number+1 of separate dense objects in the scan. Usually `N = 4`. 
        - `0` corresponds to the object with most nonzero voxels
        - `N` corresponds to the object with least nonzero voxels
        - `N` corresponds to the wooden marker
    - `color` can be one of `Red`, `Orange`, `Green` or `Blue` depending on the foam color where the panicle was located.
        - `black` corresponds to the wooden marker
    - `xcoord`, `ycoord`, `zcoord` correspond to the XYZ coordinates in the raw scan of the bottom right front corner of the bounding box containing the object of interest.

- `spikes_scanning_labels`: JPGs that depict how the colors and orientation were initially identified in each of the raw scans.

- `clean_seeds`: Resulting individual seeds from each panicle. 
    - `S001` ... `S224`: spikes corresponding to each of the 224 batches
    - `l<label>_<color>_seeds` folder for each spike within the batch
    - Seeds named `seed_<seedlabel>_<counter>_p<p>_d<d>_t<t>_o<o>_e<e>_g<g>.tif`
    - `seedlabel` goes from `0` to `N` depending on the number of seeds.
    - `counter` goes from `0` to `M`. Indicates that initially the seed was part of a larger cluster that was later broken down.
    - The rest corresponds to hyperparameters when performing the image segmentation. Refer to the corresponding Jupyter notebook in the `demeter` repository.
    
- `clean_founder_seeds`: Same information as in `clean_seeds`, except that only contains the subset of seeds corresponding to the parental generation (generation 0).
