# A Python software package for whole-brain voxel-level functional connectivity association study (Brain-wide association study)

 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Author: Weikang Gong

DPhil Student, WIN, FMRIB

Nuffield Department of Clinical Neurosciences

University of Oxford

Oxford OX3 9DU, UK

Email: weikang.gong@ndcn.ox.ac.uk

**Introduction**

Brain-wide association study is a method to analyse voxel-by-voxel resting-state fMRI data. It simply involves (1) calculating the functional connectivities between pairwise voxel across the whole brain; (2) testing the difference/correlation between functional connectivities and a phenotype of interests; (3) performing multiple comparison correction using a novel Gaussian Random Field based approach, which generalized the widely-used cluster-size inferece to functional connectivities. So intuitively, what is cluster-size inference on functional connectivities (FC), or what is FC cluster? It is just a bundle of connectivities between two voxel clusters. We are actually testing whether there are many FCs (with p-value < certain cluster-defining threshold) connecting two voxel clusters (just like in the volume analysis we test the size of the observed voxel cluster is large by chance).

This software package implement the idea of FC clusters, but also performs such a huge number of calculation (correlation matrix + GLM statistics) efficiently. It no longer has a memory requirement on your computer (but better > 16 GB), and supports both python 2.7 and 3.6, and Mac/Linux/Win system.

This is still the old matlab package. I am going to release the new python package soon in the future.









**Cite this package**
```
Gong, Weikang, et al. "Statistical testing and power analysis for brain-wide association study." Medical image analysis 47 (2018): 15-30.

Rolls, Edmund T., et al. "Functional connectivity of the anterior cingulate cortex in depression and in health." Cerebral Cortex (2018).
```


------------------------------------------------------------------------------------------------------------------------
Copyright 2018 University of Oxford

**LICENSE**

This software is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this software. If not, see http://www.gnu.org/licenses/.



**CONDITIONS OF USE**

The Software is distributed "AS IS" under the GNU General Public License (GPL, Version 3) license solely for non-commercial use. On accepting these conditions, the licensee understand that no condition is made or to be implied, nor is any warranty given or to be implied, as to the accuracy of the Software, or that it will be suitable for any particular purpose or for use under any specific conditions. Furthermore, the software authors disclaim all responsibility for the use which is made of the Software. It further disclaims any liability for the outcomes arising from using the Software.

No part of the Software may be reproduced, modified, transmitted or transferred in any form or by any means, electronic or mechanical, without the express permission of the author. The permission of the author is not required if the said reproduction, modification, transmission or transference is done without financial return, the conditions of this License are imposed upon the receiver of the product, and all original and amended source code is included in any transmitted product. You may be held legally responsible for any copyright infringement that is caused or encouraged by your failure to abide by these terms and conditions.

You are not permitted under this License to use this Software commercially. Use for which any financial return is received shall be defined as commercial use, and includes (1) integration of all or part of the source code or the Software into a product for sale or license by or on behalf of Licensee to third parties or (2) use of the Software or any derivative of it for research with the final aim of developing software products for sale or license to a third party or (3) use of the Software or any derivative of it for research with the final aim of developing non-software products for sale or license to a third party, or (4) use of the Software to provide any service to an external organisation for which payment is received. If you are interested in using the Software commercially, please contact Unitectra (http://www.unitectra.ch/en).


