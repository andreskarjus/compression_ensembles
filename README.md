**Code and data to replicate the analyses in "Compression ensembles quantify aesthetic complexity and the evolution of visual art" by Karjus et al 2023**

EPJ Data Science, 12, 21 (2023) https://doi.org/10.1140/epjds/s13688-023-00397-3
 
This code package assumes R 4.2.1, tidyverse and dplyr 1.1.0, magick 2.7.3 (an R wrapper for ImageMagick). To get started open the compression_ensembles_paper_scripts.R file in R (we recommend using RStudio), and follow the directions therein. If making use of the precomputed embeddings (stored in RData files), place them in the project folder as instructed in the script file.

The paper preprint is available here: https://arxiv.org/abs/2205.10271 
(text is updated in the EPJDS version but preprint has nicer figures, as EPJDS is unable to process high-quality images)

If you are interested in using the compression ensemble approach in general, we can also recommend looking into this optimized and more efficient Python-based version:
https://github.com/Collection-Space-Navigator/compression_ensembles

An interactive demo that includes (a subset of the) Historical (mostly wikiart) dataset is available here:
https://collection-space-navigator.github.io/
