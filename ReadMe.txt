LICENCES:
LICENCE_GNU_V3 applies to all files in the repository with the exception of the randomized data contained in dataRANDOMIZED.csv, data2016RANDOMIZED.csv and ExtendedData.pdf, which are covered by LICENCE_CC0. 

To avoid any chance that individuals in this study could be identified, the ages in the data released have been randomized by adding an integer between -5 and 5 (and then adjusted to be between 18 and 90).

To replicate the analysis on the randomized data, first run hiv.R once for each race and gender pair to generate the four MCMC chains. This step is extremely time consuming, and so shortened runs have been included by default (see comments within the file).

The files dataPlots.R and publicationPlots.R reproduce the figures.
