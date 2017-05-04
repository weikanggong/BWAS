# BWAS beta version

A Matlab-based software for brain-wide association study. To use this software, please

(1) Include all these functions in the Matlab path.

(2) Put your preprocessed fMRI data in a cell array as:

{[subject 1], [subject 2], .....[subject n]}

Each subject's data is an 2D matrix, row is the time point (t) and column is the voxels (p).

This can be done automatically by the 'BWAS_prepare.m' function. Note that 1) you should prepare a mask (the voxel-of-interests) before using this function. Voxels with non-zero elements in this mask will be used in the subsequent analysis. 2) you can save the output of this function in your hard disk.


(3) Prepare your 2D design matrix, row is the subjects, the first column is the primary variable of interests (e.g. disease status (0-1 variable), IQ (continuous variable) ), the other columns are the covariates (e.g. age, gender, motion).

This file should be prepared by yourself.

Now let's begin the analysis!

(4) After the above two steps, you will get three things: 1) the preprocessed data in a cell array, 2) the mask, 3) the design matrix with the first column as the phenotype of interest. Then, you can use the function 'BWAS_glm.m'  to perform the first step of the analysis. 

This function can 1) estimate each subjects voxel-level brain network, 2) perform connexel-wise statistical tests. The results will be automatically saved in your current matlab path as 'stat_map00*.mat'. This step may be slow.

(5) Now we can analyze the results to find significant functional connectivities and functional-connectivity clusters based on the Gaussian random field theory. This can be achieved by 'BWAS_analysis_link.m' function.



