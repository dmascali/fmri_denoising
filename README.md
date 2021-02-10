# fMRI denoising
Collection of Matlab functions for denoising fMRI data, such as resting- or steady-state data. 

- `fmri_cleaning` performs data denoising - including frequency filtering, removal of signals of no interest (e.g., motion or physiological estimates) and censoring motion-contaminated volumes – via a single linear regression model. This is the main function; the other functions allow you to create appropriate input variables for this function.

- `fmri_acompcor` extracts signals from specified regions of interest. You can either extract the average signal or extract the first n principal components (PCs;  aCompCor approach). Before extracting PCs, you can orthogonalize the data with respect to other variables that will compose the final regression model, so that the extracted PCs will be maximally predictive.  

- `fmri_rp_metrics` computes various framewise displacement metrics. 

- `fmri_censoring_mask` creates temporal masks for volume censoring (e.g., using the framewise displacement).


The code was developed for the following publication:
> **Evaluation of denoising strategies for task‐based functional connectivity: Equalizing residual motion artifacts between rest and cognitively demanding tasks.** Mascali, D., Moraschi, M., DiNuzzo, M., Tommasin, S., Fratini, M., Gili, T., Wise, R.G., Mangia, S., Macaluso, E. and Giove, F.  [*Human Brain Mapping*](https://doi.org/10.1002/hbm.25332) (2021).

