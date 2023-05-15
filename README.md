This repository contains the code accompanying our manuscript: Dissociation of vicarious and experienced rewards by coupling frequency within the same neural pathway.

The main scripts include the following: 
* `ptp_plot_original_figures.m`: Generates spike-field coherence spectra, time courses for frequency windows of interest, and time-frequency window summaries.
* `ptp_plot_raw_granger.m`: Plots time courses of spectral Granger causality for frequency windows of interest.
* `ptp_classify_sua.m`: Identifies clusters of low vs. high task-space selective cells; clusters using scores of a PCA of firing-rate metrics for each cell.

This code depends on these additional utility libraries:
* `categorical` [here](https://github.com/nfagan/categorical)
* `shared_utils` [here](https://github.com/nfagan/shared_utils)
* `bfw` [here](https://github.com/nfagan/bfw)