
# Script name: BWAS_main.py
#
# Description: Functions to run BWAS analysis using CPU
#
# Author: Weikang Gong
#
# Reference: Gong, W., Wan, L., Lu, W., Ma, L., Cheng, F., Cheng, W., Gruenewald, S. and Feng, J., 2018. Statistical testing and power analysis for brain-wide association study. Medical image analysis, 47, pp.15-30.
#
# Weikang Gong
# DPhil Student, WIN, FMRIB
# Nuffield Department of Clinical Neurosciences
# University of Oxford
# Oxford OX3 9DU, UK
# Email: weikang.gong@ndcn.ox.ac.uk or weikanggong@gmail.com
#
# Copyright 2018 University of Oxford
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#


"""
A command line tool to run BWAS analysis.

Version: 1.0.

Author: Weikang Gong.
"""


import argparse
import sys
import os


def cli_parser():
    # Create a parser. It'll print the description on the screen.
    parser = argparse.ArgumentParser(description=__doc__)
    # Add a positional argument
    parser.add_argument('-toolbox_dir', help='The absolute directory of the BWAS toolbox')
    parser.add_argument('-image_dir', help='The absolute directory of the rsfMRI images')
    parser.add_argument('-output_dir', help='The absolute directory for the BWAS outputs')
    parser.add_argument('-mask_file', help='The absolute directory of the mask file')
    parser.add_argument('-target_file', help='The absolute directory of the variable of interests (nsub * 1 numpy matrix saved in .npy format)')
    parser.add_argument('-cov_file', help='The absolute directory of the variable of interests (nsub * p numpy matrix saved in .npy format)')
    parser.add_argument('-CDT', help='Cluster defining threshold (z-value), better >= 4.5 for whole brain analysis.')
    parser.add_argument('-memory_limits', help='The maximum memory limits per CPU (GB)')
    parser.add_argument('-ncore', help='Number of CPU to use for this analysis ')
    
    return parser

parser = cli_parser()
args = parser.parse_args()

toolbox_dir=os.path.join(args.toolbox_dir,'')
output_dir=os.path.join(args.output_dir,'')
image_dir=os.path.join(args.image_dir,'')
mask_file=args.mask_file
target_file=args.target_file
cov_file=args.cov_file
CDT=float(args.CDT)
memory_limits=int(args.memory_limits)
ncore=int(args.ncore)

sys.path.append(os.path.join(os.path.abspath(toolbox_dir)))

from BWAS_cpu import BWAS_run_full_analysis


BWAS_run_full_analysis(output_dir,image_dir,mask_file,
                       toolbox_dir,target_file,cov_file,
                       CDT,memory_limits,ncore)






