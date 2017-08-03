# A Matlab software package for whole-brain voxel-level functional connectivity association study (Brain-wide association study)

The software has been tested in Matlab 2015b and high versions in the Linux CentOS 7.2 system.

A computer with 64GB memory is recommended.

**To use this software, please:**

 **1. Put all these functions in the Matlab path.**

 **2. Put the preprocessed fMRI data in a cell array as:
      {[subject 1], [subject 2], .....[subject n]}
    Each subject's data is an 2D matrix, the row is the time point (t) and column is the voxels (p).**

This can be done automatically by the 'BWAS_prepare.m' function. Note that 1) one should prepare a mask (the voxel-of-interests) before using this function. Voxels with non-zero elements in this mask will be used in subsequent analysis. Recommendation: one can use different non-zero numbers to indicate different brain regions (for example: use 1,2,...,120 for 120 regions in aal2 template). Therefore, in subsequent analysis, the locations of the significant functional connectivities will be automatically reported to you.
2) We recommend to save the output of this function in your hard disk.


**3. Prepare your 2D design matrix.** 

The row is the subjects, the first column is the primary variable of interests (e.g. disease status (0-1 variable), IQ (continuous variable) ), the other columns are the covariates (e.g. age, gender, motion). You do not need to include a     colum of 1 in the design matrix. This file should be prepared by yourself.

**Now let's begin the analysis! First, you should make a new directory for the analysis.**

**4. After the above steps, you will get three things: 1) the preprocessed data in a cell array, 2) the mask, 3) the design matrix with the first column as the phenotype of interest. Then, you can use the function 'BWAS_glm.m'  to perform the first step of the analysis.** 

This function can 1) estimate each subjects voxel-level brain network, 2) perform connexel-wise statistical tests. The results will     be automatically saved in your current matlab path as 'stat_map00*.mat'. This step may take few hours to finish.

**5. Now we can analyze the results to find significant functional connectivities and functional-connectivity clusters based on the  Gaussian random field theory. This can be achieved by 'BWAS_analysis_link.m' function. The output is a 'Link_BWAS_results.mat' file in your working directory. The descriptions within 'BWAS_analysis_link.m' illustrate the meaning of each output.**

The analysis is finished now !

**If you have any questions, please contact me (E-mail: weikanggong@gmail.com).**

