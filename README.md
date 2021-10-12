<!-- moryzae_eccdnas_manuscript_code_final -->
## Code used in manuscript describing eccDNAs in Magnaporthe oryzae

This repository provides all code used to describe eccDNAs in Magnaporthe oryzae.

The repository is organized by section of the manuscript.

Generally, each section contains:
* SLURM (.slurm) scripts that were used to first do the raw data processing on Berkeley's Savio computing cluster
* bash (.sh) and Python (.py) scripts that are called within the SLURM scripts to be used for data processing
* R (.Rmd) and iPython (.ipynb) notebooks used for final processing, statistical analyses and plotting.

The update_git.sh script as well as the two submodules are simply used by me to update and organize this github repository, please ignore them.

For acknowledgments, citations, and software versions, please refer to the manuscript:

Preprint will be posted to bioRxiv shortly

## ecc_caller

All analysis in this manuscript starts with outputs from the ecc_caller pipeline written for this manuscript. The version of this pipeline used for the manuscript is included in this repository.

A maintained version with detailed instructions for usage is posted in this separate repository:
[https://github.com/pierrj/ecc_caller](https://github.com/pierrj/ecc_caller)

<!-- LICENSE -->
## License

Code is freely available under the MIT license

<!-- CONTACT -->
## Contact

Please feel free to contact me directly with any questions, if you have any issues with the code or if you find any script that is missing.

Pierre Joubert - [@pmjoubert](https://twitter.com/pmjoubert) - pierrj@berkeley.edu

Project Link: [https://github.com/pierrj/moryzae_eccdnas_manuscript_code_final](https://github.com/pierrj/moryzae_eccdnas_manuscript_code_final)
