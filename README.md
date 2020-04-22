[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3491853.svg)](https://doi.org/10.5281/zenodo.3491853)

This repository contains all data and code necessary to run the analyses and produce the final results of the article, "Marine biodiversity and geographic distributions are independent on large scales" (Antell et al. 2020). The individual files are described below, and a visual diagram of the file sequence is supplied as the PDF named eco_release_analysis_workflow. The version numbers of important packages are listed in the Key Resources Table of the STAR Methods. Before running the code, 1) set the R working direction to the folder that contains both the /Data and /R folders, and 2) unzip the DS3 and DS4 folders so the tif and auxiliary files within are loose under /Data. The auxiliary files (.tif.xml and .tif.aux.xml) contain vital attribute information. The scripts can be run in any order, because all necessary intermediate outputs are saved in the /Data folder. The first release (v1.0.0) data and script files are archived at: https://doi.org/10.5281/zenodo.3491853. Please cite the article as follows:

Antell, G. S., Kiessling, W., Aberhan, M., & Saupe, E. E. (2020). Marine Biodiversity and Geographic Distributions Are Independent on Large Scales. Current Biology, 30(1), 115-121. https://doi.org/10.1016/j.cub.2019.10.065

/Data

- Data S1: Spreadsheet of unique spatio-temporal occurrences of brachiopod and bivalve species across
the Phanerozoic.

- Data S2: Lithochronology chart for northwest Atlantic geologic units.

- Data S3: Raster of water depth conditions for each Phanerozoic time bin.

- Data S4: Raster of substrate conditions for each Phanerozoic time bin.

- Data S5: Phanerozoic occurrence data subsampled in the assemblage framework.

- Data S6: Phanerozoic occurrence data subsampled in the taxon-tracking framework.

/R

- bin_occurrences_from_PBDB_records.R: this script downloads data from paleobiodb.org, bins the specimen records, cleans the taxonomy, converts coordinates to raster grid cells, and shortens the output spreadsheet to unique combinations of taxa, cells, and time bins. For later taxon-tracking analysis, the last part of the script also determines the environmental conditions for each occupied grid cell and the environmental preferences of each species, as described in section 5.1.1 of the Supplemental Information.

- read_data_for_bootstrap_subsampling.R: this script is sourced in both the assemblage and taxon-tracking bootstrapping scripts to generate a dataframe of time bins and format the output of bin_occurrences_from_PBDB_records.R as an object for subsampling.

ASSEMBLAGE FRAMEWORK

- bootstrap_assemblages_GSA.R: the purpose of this script is to bootstrap "assemblage-based" samples. The subsampling procedure is illustrated in Figure S1 and described in section 4.1.1 of the Supplemental Information. The script contains options to alter the subsampling parameters as in the sensitivity analyses of section 4.2. The code relies on the functions provided as function_subsample_assemblage_approach_GSA.R, function_distribution_summary_GSA.R, and function_find_seed_cells_GSA.R

- correlate_Tseries_GSA.R: this script calculates correlation coefficients between pairs of time series (e.g. geographic range size vs. species count), as described in section 4.1.2 of the Supplemental Information. The initial data are output from the "bootstrap_assemblages_GSA.R" script, and are available as "DS5_assemblage_data_subsampled.csv", which is sourced as part of the correlate_Tseries_GSA code.

TAXON-TRACKING FRAMEWORK

- categorize_enviros_from_all_PBDB_marine.R: this script takes all non-terrestrial PBDB occurrence data and uses them to characterize grid cells by environment. All taxa, at all identification levels, are included - but fossils must be "regular" not e.g. trace. The outputs of this script are 2 raster files (tif format), which are called in bin_occurrences_from_PBDB_records.R  The authors already provide these raster files, so it is unnecessary to run this script when replicating analysis. The script is provided anyway for full transparency about how the tif files were generated.

- bootstrap_across_species_GSA.R: the purpose of this script is to bootstrap species-based samples. The subsampling procedure is illustrated in Figure S2 and described in section 5.1.2 of the Supplemental Information. The script sources the function provided as function_subsample_species_approach_GSA.R

- survivor_shift_trends_GSA.R: compare mean species' range expansion/contraction during intervals of species count loss vs. increase, as described in section 5.1.4 of the Supplemental Information.

- survivor_LMM_GSA.R: The purpose of this script is to build linear mixed-effects models on species that cross stage boundaries. The initial data are output from the bootstrap_across_species_GSA.R script, and are available as DS6_species_data_subsampled.csv, which is read in here. The modeling approach is described in section 5.1.5 of the Supplemental Information.

