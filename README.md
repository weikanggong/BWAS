# BWAS: A powerful and efficient multivariate approach for voxel-level connectome-wide association studies (v0.1-beta)

## **Python code for:**

```
Weikang Gong. et al, "Statistical testing and power analysis for brain-wide association study." Medical image analysis 47 (2018): 15-30.

Edmund T. Rolls, et al. "Functional connectivity of the anterior cingulate cortex in depression and in health." Cerebral Cortex 1 (2018): 14.
```

## **Introduction:**

Brain-wide association study is a method to analyse voxel-by-voxel resting-state fMRI data. It simply involves (1) calculating the functional connectivities between pairwise voxel across the whole brain; (2) testing the difference/correlation between functional connectivities and a phenotype of interests; (3) performing multiple comparison correction using a novel Gaussian Random Field-based approach, which generalized the widely-used cluster-size inferece to functional connectivities. So intuitively, what is cluster-size inference on functional connectivities (FC), or what is FC cluster? It is just a bundle of connectivities between two voxel clusters. We are actually testing whether there are many FCs (with p-value < certain cluster-defining threshold) connecting two voxel clusters (just like in the volume analysis we test the size of the observed voxel cluster is large by chance). This software package implement the idea of FC clusters, but also performs such a huge number of calculations (correlation matrix + GLM statistics) efficiently. It no longer has a memory requirement on your computer (but better $>$ 16 GB), and supports both python 2.7 and 3.6, and Linux/Mac/Windows system.

This is the development version of BWAS, bug report is wellcome!

## **Requirement:**
1. System: Linux/Mac/Windows;
2. Python 2.7 or 3.6 (Anaconda is recommended);
3. Python modules: copy, glob, numpy, scipy, matplotlib, nilearn, nibabel, joblib, multiprocessing, PyPDF2;


## **Data structure and required files:**
1. Toolbox directory: The absolute directory of the BWAS code;
2. fMRI data: Please put all your rfMRI data in a directory. The software will read data in alphabet order.
3. variable of interest file: One column. The file format should be either a ".txt" file or a ".npy" file, with each row representing a subject and column representing a variable.
4. covariates file: Multiple columns. The file format should be either a ".txt" file or a ".npy" file, with each row representing a subject and each column representing a variable.
5. mask_file: a binary mask (.nii.gz or .nii format) of your fMRI data;
6. CDT: cluster defining threshold (z-value), usually >=5;
7. Memory limits: the maxmum number of memory to use per CPU (in GB);
8. Number of cores: usually the more the faster; Memory limits * Number of cores must be < your total avaliable memory.
9. Output directory: the absolute directory to save all the outputs.

## **How to use this package:**
1. All the source code is in the file: **BWAS_cpu.py**
2. To run it in command line, please use the file: BWAS_main.py; You can type: **python BWAS_main.py -h** to see the help;
3. To run it in GUI, please use the file: BWAS_gui.py; You can type: **python BWAS_gui.py** to open the gui, the input should be the same as BWAS_main.py. After enter all the things, press "run BWAS interactively" or "run BWAS in background" to perform the analysis (Do not forget to press the save button after you enter CDT,Memory limits and Number of cores).


## **Outputs:**

In the Output directory, 
1. BWAS_result_CDT=?.mat: the infomation of all FC with |Z|> CDT, including the matrix coordinates and MNI coordinates and the corresponding Z-statistics of each FC.
2. 'MA_CDT=?_x.jpg','MA_CDT=?_y.jpg','MA_CDT=?_z.jpg': the figures of the number of significant FCs in each voxel;
3. 'MA_CDT=??.nii.gz': the number of significant FCs in each voxel;
4.  BWAS_FC_clusters_plot.pdf: Each figure represent a FC cluster, connecting the left voxel-cluster and the right voxel-cluster, the values on the voxels shows the number of FCs connecting it that pass FC cluster based correction;
5. ./FC_cluster_files/FC_cluster_?+'_1.nii.gz, ./FC_cluster_files/FC_cluster_?+'_2.nii.gz: the nifti files to generate the above pictures.
6. brain_gl_negative.txt and brain_gl_positive.txt: the file you can put into braingl software (https://code.google.com/archive/p/braingl/) to visualize the voxel-by-voxel FCs. You should use "braingl_bg.nii.gz" as background. 


## **Question or report bug:**

Author: Weikang Gong (FMRIB Analysis group, NDCN, WIN, University of Oxford)

Email: weikang.gong@ndcn.ox.ac.uk; weikanggong@gmail.com


You may also find sKPCR is useful: https://github.com/weikanggong/sKPCR.

