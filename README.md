# A Matlab software package for whole-brain voxel-level functional connectivity association study (Brain-wide association study)

The software has been tested in Matlab 2015b and higher versions in the Linux CentOS 7.2 system.

A computer with 64GB memory is recommended.

**To use this software, please prepare the following 4 things:**

 1. Put all these functions in the Matlab path.

 2. Put all the fMRI images in a folder. Each image must be in .nii.gz format.
 
 3. Generate a design matrix in *.txt format. The number of rows equal to the number of subjects in the analysis, the order is the same as the images' order when typing 'dir('*.nii.gz')' in Matlab. The first column of the design matrix is the phenotype of interest (e.g. disease status (0-1 variable), IQ (continuous variable)), and other columns are the covariates (e.g. age, gender, motion).
 
 4. Generate a brain mask image in *.nii or *.nii.gz format. The mask is with the same size as the fMRI data. Non-zeros elements in this mask indicates the voxels one wants to analyse.
 
 **Then you can simply use the function**
 ```
 BWAS_main(result_dir,image_dir,design_dir,mask_dir,CDT,FWER_p)
 ```
 **to perform all the analysis. If the analysis is interupted, simply rerun it. It will automatically continue the steps that have not been performed.**
 
 
 **The inputs of the 'BWAS_main' function are quite simple:**
 ```
result_dir: it is the absolute directory of a folder to save the outputs of the analysis.
image_dir: it is the absolute directory of a folder of your images.
design_dir: it is the absolute directory of the design matrix in *.txt format.
mask_dir: it is the absolute directory of the mask file in *.nii or *.nii.gz format.
CDT: Cluster-definding threshold of functional connectivity clusters, default is Z=5 if one performs a whole-brain analysis.
FWER_p: the p-value threshold of peak-level inference, default is 0.05.
 ```
 **The outputs of the 'BWAS_main' function are illustrated in the 'BWAS_analysis_link' function. The most important things are:**
 
 1. The 6D coordinates of the functional connectivities clusters exceeding the CDT, and the corresponding FWER and FDR corrected p-values of their size. See 'Link_BWAS_results.mat' for details.
 
 2. The number of significant functional connectivities connecting each voxel. See 'peak_MA.nii.gz' and 'cluster_MA.nii.gz' for details.
 
**An example of design matrix (4 subjects, 1 phenotype of interest and 3 covariates):**
```
0 1.5 2 3
0 2.2 5 4
1 2   3 1
1 2   4 5
```


**If you have any questions, please contact me (E-mail: weikanggong@gmail.com).**

