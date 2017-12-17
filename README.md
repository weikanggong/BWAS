# A Matlab software package for whole-brain voxel-level functional connectivity association study (Brain-wide association study)

 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Copyright (C) 2017 Weikang Gong

The content in this package is licensed under a Creative Commons Attribution-Noncommercial-Share Alike 3.0 United States License.

**You are free to:**

Share — copy and redistribute the material in any medium or format.

Adapt — remix, transform, and build upon the material.

**Under the following terms:**

Attribution — You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.

NonCommercial — You may not use the material for commercial purposes.

ShareAlike — If you remix, transform, or build upon the material, you must distribute your contributions under the same license as the original.

No additional restrictions — You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

**Notices:**

You do not have to comply with the license for elements of the material in the public domain or where your use is permitted by an applicable exception or limitation.

No warranties are given. The license may not give you all of the permissions necessary for your intended use. For example, other rights such as publicity, privacy, or moral rights may limit how you use the material.
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


The software has been tested in Matlab 2015b and higher versions in the Linux CentOS 7.2 system.

A computer with 64GB memory or more is recommended.

**To use this software, please prepare the following 4 things:**

 1. Put all these functions in the Matlab path.

 2. Put all the fMRI images in a folder. Each image must be in .nii.gz format.
 
 3. Generate a design matrix in *.txt format. The number of rows equal to the number of subjects in the analysis, the order is the same as the images' order when typing 'dir('*.nii.gz')' in Matlab. The first column of the design matrix is the phenotype of interest (e.g. disease status (0-1 variable), IQ (continuous variable)), and other columns are the covariates (e.g. age, gender, motion).
 
 4. Generate a brain mask image in *.nii or *.nii.gz format. The mask is with the same size as one volume of the fMRI data. Non-zeros elements in this mask indicates the voxels one wants to analyse. **Important: make sure that all the fMRI time series within the mask are not constant value.**
 
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

