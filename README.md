# widefield-cerebellar-calcium-imaging
code used in https://www.biorxiv.org/content/10.1101/2023.10.05.561090v1



These analysis scripts were written using Matlab 2019b, and has not been tested on any other versions.

Download and add scripts to current directory to run.

Summary of scripts and expected output:

1) preprocessing_wheel.m concatenates calcium imaging stacks during spontaneous locomotion paradigm, performs preprocessing, and compresses the resulting preprocessed imaging using singular value decomposition
2) preprocessing_reaching_MSI.m same as above, but for imaging data obtained during the cued reaching paradigm
3) run_ICA_on_compressed.m runs sICA on resulting compressed data files
4) ICA_manual_threshold.m script for IC thresholding, allowing for the manual discard of ICs corresponding to artifacts e.g. vasculature
5) IC_dff.m takes all ICs generated from the previous script, extracts df/f for each IC domain, and saves the results
6) locomotion_segment_collectdff.m computes average df/f for each IC for each locomotion epoch
7) locomotion_epochs_corrmaps_SOM.m constructs spatial correlation maps for each individual mouse, each individual recording day, for each locmotion epoch within that recording day, all for somatic ICs only
8) locomotion_epochs_corrmaps_DEND.m same as above but for dendritic ICs
9) bin_av_corrmaps_locomotion.m averages and bins the invidual spatial correlation maps computed in the previous script to generate population maps
10) reaching_epochs_corrmaps_SOM.m generates spatial correlation maps for individual mice, individual recording days, during the reaching paradigme pochs, all for somatic ICs only
11) reaching_epochs_corrmaps_DEND.m same as above but for dendritic ICs
12) bin_av_corrmaps_reaching.m averages and bins the invidual spatial correlation maps computed in the previous script to generate population maps
13) reaching_catch_epochs_corrmaps_SOM.m generates spatial correlation maps for individual mice, individual recording days, during the catch paradigm epochs, all for somatic ICs only
14) reaching_catch_epochs_corrmaps_DEND.m same as bove, but for dendritic ICs
15) bin_av_corrmaps_catch.m averages and bins the invidual spatial correlation maps computed in the previous script to generate population maps
16) computezscore_standardreach_truereach.m goes through each IC for each mouse and computes the z score of modulation for successful reaches, thresholds those ICs, and constructs day level maps of thresholded IC modulation
17) computezscore_standardreach_noreach.m same as above, but for thresholded IC modulation during failed reaches
18) computezscore_catchtrial.m same as above, but for the catch trial
19) plot_collapsed_z_maps_reaching.m takes the maps constructed from the three previous scripts and averages them to create time series population maps for each reaching period of interest
20) sliding_window_reachnoreach_som.m performs sliding window correlation analysis for somatic ICs, aligned to the auditory tone, when the mouse either reached or failed to reach during the standard reaching paradigm
21) sliding_window_reachnoreach_dend.m same as above, but for dendritic ICs
22) sliding_window_catchreach_som.m performs sliding window correlation analysis for somatic ICs, aligned to the auditory tone, when the mouse either reached or during the catch trial
23) sliding_window_catchreach_dend.m same as above, but for dendritic ICs

