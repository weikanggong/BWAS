
# Script name: BWAS_gui.py
#
# Description: Functions to run BWAS analysis using CPU in GUI
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




import warnings
warnings.filterwarnings("ignore")
import sys
import os

if sys.version_info[0]==2:
    from Tkinter import *
    from ttk import *
    from tkFileDialog import *
else:
    from tkinter import *
    from tkinter.ttk import *
    from tkinter.filedialog import *




def def_toolbox_dir():
    global toolbox_dir
    toolbox_dir = askdirectory()
    entry_toolbox.delete(0, END)
    entry_toolbox.insert(0, toolbox_dir)
    return toolbox_dir

def def_data_folder():
    global data_folder
    data_folder=askdirectory()
    entry_data.delete(0, END)
    entry_data.insert(0, data_folder)
    return data_folder

def def_out_dir():
    global out_dir
    out_dir = askdirectory()
    entry_out.delete(0, END)
    entry_out.insert(0, out_dir)
    return out_dir

def def_behavior_data_file1():
    global behavior_data1
    behavior_data1 = askopenfilename(filetypes=[("NPY or text files","*.npy *.txt")])

    entry_b1.delete(0, END)
    entry_b1.insert(0, behavior_data1)
    return behavior_data1

def def_behavior_data_file2():
    global behavior_data2
    behavior_data2 = askopenfilename(filetypes=[("NPY or text files","*.npy *.txt")])
    
    entry_b2.delete(0, END)
    entry_b2.insert(0, behavior_data2)
    return behavior_data2

def def_image_data_file():
    global image_data
    image_data = askopenfilename(filetypes=[("Nifti files","*.gz *.nii")])
    
    entry_mask.delete(0, END)
    entry_mask.insert(0, image_data)
    return image_data

global CDT
CDT=str(5)
def callback3():
    global CDT
    CDT = entry3.get()
    entry3.delete(0,END)
    entry3.insert(0,CDT)
    return CDT

global memory_limit
memory_limit=str(8)
def callback4():
    global memory_limit
    memory_limit = entry4.get()
    memory_limit = memory_limit
    entry4.delete(0,END)
    entry4.insert(0,memory_limit)
    return memory_limit

global ncore
ncore=str(1)
def callback5():
    global ncore
    ncore = entry5.get()
    ncore = ncore
    entry5.delete(0,END)
    entry5.insert(0,ncore)
    return ncore


def gui_calls_BWAS():
    #import pdb;pdb.set_trace()
    os.system('chmod +x '+ os.path.join(toolbox_dir,'BWAS_main.py') +
              '&& python '+os.path.join(toolbox_dir,'BWAS_main.py') +
              ' -toolbox_dir ' + toolbox_dir +
              ' -image_dir ' + data_folder + 
              ' -output_dir ' + out_dir + 
              ' -mask_file ' + image_data + 
              ' -target_file ' + behavior_data1 + 
              ' -cov_file ' + behavior_data2 +
              ' -CDT ' + CDT +
              ' -memory_limits ' + memory_limit +
              ' -ncore ' + ncore )
    return 

def gui_calls_BWAS_back():
    #import pdb;pdb.set_trace()
    os.system('chmod +x '+ os.path.join(toolbox_dir,'BWAS_main.py') +
              '&& nohup python '+os.path.join(toolbox_dir,'BWAS_main.py') +
              ' -toolbox_dir ' + toolbox_dir +
              ' -image_dir ' + data_folder + 
              ' -output_dir ' + out_dir + 
              ' -mask_file ' + image_data + 
              ' -target_file ' + behavior_data1 + 
              ' -cov_file ' + behavior_data2 +
              ' -CDT ' + CDT +
              ' -memory_limits ' + memory_limit +
              ' -ncore ' + ncore + ' &')
    return 


root = Tk() # create a top-level window
 
master = Frame(root, name='master') # create Frame in "root"
master.pack(fill=BOTH) # fill both sides of the parent
 
root.title('BWAS - voxel-level functional connectivity analysis') # title for top-level window
# quit if the window is deleted
root.protocol("WM_DELETE_WINDOW", master.quit)
 
nb = Notebook(master, name='nb') # create Notebook in "master"
nb.pack(fill=BOTH, padx=2, pady=3) # fill "master" but pad sides
 
#-->INPUT DATA TAB
f1 = Frame(nb, width=600, height=250)
f1.pack(fill=X)
nb.add(f1, text="Set up and run BWAS analysis") # add tab to Notebook


Label(f1,text="Select Toolbox directory").grid(row=0, column=0, sticky='e')
entry_toolbox = Entry(f1, width=50)#, textvariable=folder_path)
entry_toolbox.grid(row=0,column=1,padx=2,pady=2,sticky='we',columnspan=25)
Button(f1, text="Browse", command=def_toolbox_dir).grid(row=0, column=27, sticky='ew', padx=8, pady=4)

Label(f1,text="Select input rsfMRI data directory").grid(row=1, column=0, sticky='e')
entry_data = Entry(f1, width=50)#, textvariable=folder_path)
entry_data.grid(row=1,column=1,padx=2,pady=2,sticky='we',columnspan=25)
Button(f1, text="Browse", command=def_data_folder).grid(row=1, column=27, sticky='ew', padx=8, pady=4)

Label(f1,text="Select output directory").grid(row=2, column=0, sticky='e')
entry_out = Entry(f1, width=50)#, textvariable=folder_path)
entry_out.grid(row=2,column=1,padx=2,pady=2,sticky='we',columnspan=25)
Button(f1, text="Browse", command=def_out_dir).grid(row=2, column=27, sticky='ew', padx=8, pady=4)

Label(f1,text="Select mask file (.nii.gz or .nii)").grid(row=3, column=0, sticky='e')
entry_mask = Entry(f1, width=50)#, textvariable=folder_path)
entry_mask.grid(row=3,column=1,padx=2,pady=2,sticky='we',columnspan=25)
Button(f1, text="Browse", command=def_image_data_file).grid(row=3, column=27, sticky='ew', padx=8, pady=4)

Label(f1,text="Select variable of interest file (.txt or .npy)").grid(row=4, column=0, sticky='e')
entry_b1 = Entry(f1, width=50)#, textvariable=folder_path)
entry_b1.grid(row=4,column=1,padx=2,pady=2,sticky='we',columnspan=25)
Button(f1, text="Browse", command=def_behavior_data_file1).grid(row=4, column=27, sticky='ew', padx=8, pady=4)

Label(f1,text="Select covariates file (.txt or .npy)").grid(row=5, column=0, sticky='e')
entry_b2 = Entry(f1, width=50)#, textvariable=folder_path)
entry_b2.grid(row=5,column=1,padx=2,pady=2,sticky='we',columnspan=25)
Button(f1, text="Browse", command=def_behavior_data_file2).grid(row=5, column=27, sticky='ew', padx=8, pady=4)

Label(f1,text="Cluster Defining Threshold (z-value, Default = 5)").grid(row=6, column=0, sticky='e')
entry3 = Entry(f1, width=10)#, textvariable=dum_str)
entry3.grid(row=6,column=1,padx=2,pady=2,sticky='we',columnspan=10)
Button(f1, text="Save", command=callback3).grid(row=6, column=27, sticky='ew', padx=8, pady=4)


Label(f1,text="Memory limit (per cpu in GB, Default = 8)").grid(row=7, column=0, sticky='e')
entry4 = Entry(f1, width=10)#, textvariable=dum_str)
entry4.grid(row=7,column=1,padx=2,pady=2,sticky='we',columnspan=10)
Button(f1, text="Save", command=callback4).grid(row=7, column=27, sticky='ew', padx=8, pady=4)

Label(f1,text="Number of CPU to use (Default = 1)").grid(row=8, column=0, sticky='e')
entry5 = Entry(f1, width=10)#, textvariable=dum_str)
entry5.grid(row=8,column=1,padx=2,pady=2,sticky='we',columnspan=10)
Button(f1, text="Save", command=callback5).grid(row=8, column=27, sticky='ew', padx=8, pady=4)


Button(f1, text="Run BWAS interactively", command=gui_calls_BWAS).grid(row=11, column=10, sticky='ew', padx=8, pady=4)

Button(f1, text="Run BWAS in background", command=gui_calls_BWAS_back).grid(row=14, column=10, sticky='ew', padx=8, pady=4)

Button(f1, text="KILL GUI", command=master.quit).grid(row=17, column=10, sticky='ew', padx=8, pady=4)


 
# start the app
if __name__ == "__main__":
    master.mainloop() # call master's Frame.mainloop() method.
    #root.destroy() # if mainloop quits, destroy window

