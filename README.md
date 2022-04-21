# Readme

This repository provides data and code to reproduce the figures of **"Contrasting state-dependent effects of natural forcing on global and local climate variability"** (Ellerhoff et al., 2022) accepted in *Geophysical Research Letters*. The presented code is released under a Zenodo identifier.

**Authors:** Beatrice Ellerhoff, Moritz J. Kirschner, Elisa Ziegler, Kira Rehfeld

**Responsibility for this repository:** Beatrice Ellerhoff ([@bellerhoff](https://github.com/bellerhoff)) and Kira Rehfeld ([@krehfeld](https://github.com/krehfeld)). The code builds on previous analysis by Moritz J. Kirschner ([@cellador](https://github.com/cellador)). Elisa Ziegler ([@elisaziegler](https://github.com/elisaziegler)) helped to conduct the TransEBM simulations and their analysis (https://github.com/paleovar/TransEBM). 

**Related persons**: Max D. Holloway, Louise Sime

Please see `./license.md` for terms of use. This repository contains the **maintained code and processed data** to create the figures of *Ellerhoff et al. (2022)*. The model simulation data can be found on [Zenodo](https://doi.org/10.5281/zenodo.6074747) under the identifier 10.5281/zenodo.6074747 . 

## Organisation of this repository

The repository contains scripts (e.g., `F2.R`), input data (`./data`), and additional information (e.g., `license.md`). It is organized in a folder structure, where each folder comprises a specific analysis. To reproduce a figure from *Ellerhoff et al. (2022)* run the script that is named as the figure (Figure `X` of the main manuscript is denoted by `FX.R` and Figure `X` of the supplementary materials is denoted by `FSX.R`). The main directory contains  scripts (`.R`-files) and files that define useful functions and plotting parameters.

scripts | description
---- | ----------
`F1.R`- `F4.R` | Scripts to reproduce the figures of the **main manuscript**.
`FS1.R`- `FS11.R`| Scripts to reproduce the figures of the **supplemental information**.

directories | description
---- | ----------
`./01_characteristics/` | Contains the analysis for characterizing the input data, i.e. model runs.  
`./02_surface_climate_response/` | Contains all analysis of the surface climate response to volcanism.
`./03_spectral_analysis/` | Contains the spectral analysis of the effects from natural forcing on global and regional variability. 
`./04_variance_ratios/` | Contains the variance ratio analysis for comparison of simulated to observed variance.
`./data/` | Contains pre-processed input data useful for plotting. 
`./data/HadCM3/` | Empty. Please load the raw_data from the [Zenodo](https://doi.org/10.5281/zenodo.6074747) identifier 10.5281/zenodo.6074747 

additional files | description
---- | ----------
`.gitignore` | Information for GIT version control to not add several file extensions to version control (e.g. `*.png`, `*.pdf`)
`license.md`/ `license.html` | Licensing information
`readme.md` | General README

## Prerequisites

Our code requires the following [R](https://www.r-project.org/) packages:

- `stats`
- `utils`
- `polyclip` (might additionally require [`Rtools`](https://cran.r-project.org/))
- `ncdf4`
- `RColorBrewer`
- `maps` (can be obtained from https://github.com/krehfeld/nest, using `devtools::install_github()`)
- `PaleoSpec` (can be obtained from https://github.com/EarthSystemDiagnostics/paleospec, using `devtools::install_github()`)

Creating some of the figures additionally requires (separately loaded in the corresponding `.R` scripts):

- `dplyr`
- `openxlsx`

Generating the input data once the raw data is available on Zenodo (after publication) additionally requires:

- `ncdf4`
- `boot`
- `parallel`
- `graphics`

## Data references

We thank the research groups for producing and making available their data from paleoclimate, observations, and forcing reconstructions.

Model simulations were carried out using:
- the model TransEBM **E. Ziegler and K. Rehfeld**, described in *TransEBM v. 1.0: description, tuning, and validation of a transient model of the Earth's energy balance in two dimensions*, Geoscientific Model Development (2021). Please find the current version of the model at https://github.com/paleovar/TransEBM. The first version is licensed under the GNU General Public License Version 3 and archived on Zenodo (https://doi.org/10.5281/zenodo.3941311, Ziegler and Rehfeld, 2020).
- version 3 of the Hadley Center Coupled Model, HadCM3, described in **P. Valdes et al.**, *The BRIDGE HadCM3 family of climate models: HadCM3@Bristol v1.0*, Geoscientific Model Development (2017) and *J. Tindall, L. Sime, P. Valdes* Journal of Geophysical Research (Atmospheres) (2009) https://doi.org/10.1029/2008JD010825 . The model simulation data is published on [Zenodo](https://doi.org/10.5281/zenodo.6074747) under the identifier 10.5281/zenodo.6074747 
- the Archer UK National Supercomputing Services (www.archer.ac.uk) for HadCM3.

Paleoclimate, observational and reanalysis data were obtained from:

- **Rehfeld et al.**, *Global patterns of declining temperature variability from the Last Glacial Maximum to the Holocene*, Nature (2018) https://doi.org/10.1038/nature25454
- **Rayner et al.**. *Global analyses of sea surface temperature, sea ice, and night marine air temperature since the late nineteenth century*, Journal of Geophysical Research D: Atmospheres (2003) (downloaded latest version of HadISST dataset 11/2019)
- **PAGES 2k Consortium**, "A global multiproxy database for temperature reconstructions of the Common Era." Scientific data 4 (2017).https://dx.doi.org/10.1038%2Fsdata.2017.88

Forcing reconstructions were obtained from:
- **G. A. Schmidt et al.**, *Climate forcing reconstructions for use in PMIP simulations of the Last Millennium (v1.1)*, Geoscientific Model Development (2012)
- **T. J. Crowley and M. B. Unterman**, *Technical details concerning development of a 1200 yr proxy index for global volcanism*, Earth System Science Data (2013)
- **F. Steinhilber et al.**, *Total solar irradiance during the Holocene*, Geophysical Research Letters (2009)
- **Y. Wang et al.**, *Modeling the Sun’s Magnetic Field and Irradiance since 1713*, The Astrophysical Journal (2005)

---
We acknowledge the [R Core team](https://www.R-project.org/) and all package developers of packages used in this study. We thank them for their time and dedication to provide R and the packages to the public. Please see `citation()` for details on the R Core Team and `citation("pkgname")` for details on the developers of individual packages.

The study *Ellerhoff et al. (2022)* has been supported by funds of the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation), [Project No. 395588486](https://gepris.dfg.de/gepris/projekt/395588486?context=projekt&task=showDetail&id=395588486&) and [Project No. 316076679](https://gepris.dfg.de/gepris/projekt/316076679?context=projekt&task=showDetail&id=316076679&), by the [PalMod](https://www.palmod.de/) project (subproject no. 01LP1926C), and the [Heinrich-Böll-Stiftung](https://boell.de/). HadCM3 simulations were performed on the Archer UK National Supercomputing Services. The study benefited from discussions within the [CVAS](https://pastglobalchanges.org/science/wg/cvas/intro) working group, a working group of the [Past Global Changes (PAGES)](https://pastglobalchanges.org/pal) project. We thank the members of the [Earth's climate and environmental dynamics](https://www.iup.uni-heidelberg.de/en/research/paleoclimate-dynamics) and the [SPACY](https://uni-tuebingen.de/climatology/) groups for discussion at different stages of the manuscript preparation. 

Please report bugs to the authors (beatrice.ellerhoff(at)stud.uni-heidelberg.de, kira.rehfeld(at)uni-tuebingen.de).

*Beatrice Ellerhoff and Kira Rehfeld, April 2022*
