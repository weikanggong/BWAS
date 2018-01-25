# A Matlab software package for whole-brain voxel-level functional connectivity association study (Brain-wide association study)

 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Copyright (C) 2017 Weikang Gong

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

**Cite this package**
```
Gong, Weikang, et al. "Statistical testing and power analysis for brain-wide association study." bioRxiv (2017): 089870.
```

**If you have any questions, please contact me (E-mail: weikanggong@gmail.com).**

**LICENSE**

This software is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this software. If not, see http://www.gnu.org/licenses/.



**CONDITIONS OF USE**

The Software is distributed "AS IS" under the GNU General Public License (GPL, Version 3) license solely for non-commercial use. On accepting these conditions, the licensee understand that no condition is made or to be implied, nor is any warranty given or to be implied, as to the accuracy of the Software, or that it will be suitable for any particular purpose or for use under any specific conditions. Furthermore, the software authors disclaim all responsibility for the use which is made of the Software. It further disclaims any liability for the outcomes arising from using the Software.

No part of the Software may be reproduced, modified, transmitted or transferred in any form or by any means, electronic or mechanical, without the express permission of the author. The permission of the author is not required if the said reproduction, modification, transmission or transference is done without financial return, the conditions of this License are imposed upon the receiver of the product, and all original and amended source code is included in any transmitted product. You may be held legally responsible for any copyright infringement that is caused or encouraged by your failure to abide by these terms and conditions.

You are not permitted under this License to use this Software commercially. Use for which any financial return is received shall be defined as commercial use, and includes (1) integration of all or part of the source code or the Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use of the Software or any derivative of it for research with the final aim of developing software products for sale or license to a third party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software products for sale or license to a third party, or (4) use of the Software to provide any service to an external organisation for which payment is received. If you are interested in using the Software commercially, please contact Unitectra (http://www.unitectra.ch/en).


