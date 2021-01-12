# fMRI denoising
Collection of Matlab functions for denoising fMRI data, such as resting- or steady-state data. 

- `fmri_cleaning` performs data denoising - including frequency filtering, removal of signals of no interest (e.g., motion or physiological estimates) and censoring motion-contaminated volumes – via a single linear regression model. This is the main function; the other functions allow you to create appropriate input variables for this function.

- `fmri_compcor` extracts signals from specified regions of interest. You can extract the average signal, the median signal or the principal components (PCs;  *aCompCor* or *tCompcor*). Before extracting PCs, you can orthogonalize the data with respect to other variables that will compose the final regression model, so that the extracted PCs will be maximally predictive.  

- `fmri_rp_metrics` computes various framewise displacement metrics. 

- `fmri_censoring_mask` creates temporal masks for volume censoring (e.g., using the framewise displacement).

