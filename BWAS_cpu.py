
# Script name: BWAS_cpu.py
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




import time
import os 
from copy import deepcopy
import glob
import sys

import numpy as np

from scipy.spatial.distance import cdist
from scipy import sparse
from scipy.stats import norm
import scipy
import scipy.io as sio
import scipy.signal
from scipy import misc

import matplotlib as mpl
mpl.use('Agg')
from nilearn import plotting

import nibabel as nib

from joblib import Parallel, delayed
import multiprocessing
from PyPDF2 import PdfFileMerger


def BWAS_correlation(fMRI_2D_1,fMRI_2D_2):
    # fMRI_2D_1 and fMRI_2D_2 are both time * voxel (t * p1 and t * p2) matrices
    # This function return the fisher z transformed correlation matrix (p1 * p2)
    
    fMRI_2D_1 = (fMRI_2D_1 - fMRI_2D_1.mean(axis=0)) / fMRI_2D_1.std(axis=0)
    fMRI_2D_2 = (fMRI_2D_2 - fMRI_2D_2.mean(axis=0)) / fMRI_2D_2.std(axis=0)
    r=np.dot(np.transpose(fMRI_2D_1),fMRI_2D_2) / fMRI_2D_1.shape[0]
    
    return r

def BWAS_fisher_z(r):
    r=0.5 * np.log(np.divide(1+r,1-r))
    return r


def BWAS_regression_online1(X,Y,st,en,Betas0):
    
    ss=np.linalg.inv(np.dot(X.T,X))
    Betas1=Betas0+np.dot(np.dot(ss,X[st:en,:].T),Y)
    
    return Betas1

def BWAS_regression_online2(X,Y,st,en,Beta_glm,Sigma0):

    #caculate current residual
    Res=Y-np.dot(X[st:en,:],Beta_glm)

    #update sigma
    Sigma_glm=Sigma0+np.sum(np.square(Res),axis=0)
    return Sigma_glm


def BWAS_regression_online3(X,Beta_glm,Sigma_glm):
    
    #contrast
    contrast=np.hstack((np.ones((1,1),dtype='float32'),np.zeros((1,X.shape[1]-1),dtype='float32'))).T
    
    ss=np.linalg.inv(np.dot(X.T,X))

    #df
    df=X.shape[0]-X.shape[1]
    
    #tstat
    Tstat=np.divide(np.dot(Beta_glm.T,contrast),np.dot(np.sqrt(Sigma_glm/df).T,np.sqrt(np.dot(np.dot(contrast.T,ss),contrast))))

    return Tstat

def BWAS_est_fwhm(fMRI):
    #4d fmri + 3d mask
    #output fwhm of each time point
    #mask=nib.load('/Users/wgong/Documents/MATLAB/hcp_bwas/hcp_bwas_mask_final_4mm.nii.gz').get_data()
    #fMRI=nib.load('/Users/wgong/Documents/MATLAB/hcp_bwas/100206_rest1.nii.gz').get_data()
    
    fMRI=scipy.signal.detrend(fMRI,3)
    
    if len(fMRI.shape)==3:
        [n1,n2,n3]=fMRI.shape
        n4=1
    else:   
        [n1,n2,n3,n4]=fMRI.shape
        
    fMRI=fMRI*fMRI/fMRI
    
    dx = np.diff(fMRI,1,0); 
    varX=np.nanvar(np.reshape(dx,[(n1-1)*n2*n3,n4]),axis=0)

    dy = np.diff(fMRI,1,1); 
    varY=np.nanvar(np.reshape(dy,[n1*(n2-1)*n3,n4]),axis=0)

    dz = np.diff(fMRI,1,2); 
    varZ=np.nanvar(np.reshape(dz,[n1*n2*(n3-1),n4]),axis=0)

    varXYZ=(varX*varY*varZ)**(1/3.0);
        
    varImg=np.zeros((n4,))
    for i in range(0,n4):
        tmp=fMRI[:,:,:,i]
        varImg[i]=np.nanvar(tmp)

    fwhm=np.real(np.sqrt(-2*np.log(2)/np.log(1-varXYZ/2/varImg)))
		    
    return fwhm

def BWAS_prepare(imgs_abs_dir,mask,nargout):
#        BWAS_prepare(imgs_abs_dir,mask,name)
#Input   imgs_abs_dir: list, one cell one absolute path of preprocessed fmri image
#        mask: 3D matrix, non-zero elements are the voxel you want to use
#        in the analysis.
#Output  images: the data can be used in subsequent analysis.


    tmp1= '.nii.gz' in imgs_abs_dir[0][0]
    #tmp2= '.dtseries.nii' in imgs_abs_dir[0][0]

    nrun=len(imgs_abs_dir[0])
    nsub=len(imgs_abs_dir)
    mask1=np.reshape(mask,[mask.shape[0]* mask.shape[1]* mask.shape[2]])
    
    list_of_arrays=[np.array(a) for a in range (0,nsub)]
    images = deepcopy(list_of_arrays)
    
    
    fwhm=np.zeros((nsub,nrun))
    
    if tmp1:
        for i in range(0,nsub):

            image1=None
            for j in range(0,nrun):
                print('Reading image '+imgs_abs_dir[i][j]+'...')
                #load data
                img = nib.load(imgs_abs_dir[i][j])        
                data = np.float32(img.get_data())
                #estimate fwhm
                if nargout==2:
                    print('Estimating Smoothness...')
                    fwhm[i,j]=np.mean(BWAS_est_fwhm(data))
                    print('FWHM = '+str(fwhm[i,j])+' voxels')

                #to 2d matrix
                img2d=np.reshape(data,[data.shape[0]* data.shape[1]* data.shape[2],data.shape[3]])
                img2d=np.transpose(img2d[mask1!=0,:])
                #connect by time
                
                if j==0:
                    image1=(img2d - img2d.mean(axis=0)) / img2d.std(axis=0)
                else:    
                    image1=np.vstack((image1,(img2d - img2d.mean(axis=0)) / img2d.std(axis=0)))
                
                
                print('Done...')
           
            images[i]=image1
                   
        if nargout==2:
            fwhm=np.mean(fwhm)
            return image1,fwhm 
        if nargout==1:
            return image1

def BWAS_prepare_parallel(imgs_abs_dir,mask,ncore,nargout):
#        BWAS_prepare(imgs_abs_dir,mask,name)
#Input   imgs_abs_dir: list, one cell one absolute path of preprocessed fmri image
#        mask: 3D matrix, non-zero elements are the voxel you want to use
#        in the analysis.
#Output  images: the data can be used in subsequent analysis.


    
    def readImage(img_name,mask,nrun,nargout):
        
        fwhm=np.zeros((nrun,))
        for j in range(0,nrun):
            print('Reading image '+img_name[j]+'...')
            #load data
            img = nib.load(img_name[j])        
            data = np.float32(img.get_data())
            #estimate fwhm
            if nargout==2:
                print('Estimating Smoothness...')
                fwhm[j]=np.mean(BWAS_est_fwhm(data))
                print('FWHM = '+str(fwhm[j])+' voxels')

            #to 2d matrix
            img2d=np.reshape(data,[data.shape[0]* data.shape[1]* data.shape[2],data.shape[3]])
            img2d=np.transpose(img2d[mask1!=0,:])
            #connect by time
            
            if j==0:
                image1=(img2d - img2d.mean(axis=0)) / img2d.std(axis=0)
            else:    
                image1=np.vstack((image1,(img2d - img2d.mean(axis=0)) / img2d.std(axis=0)))
            
            
            print('Done...')
        if nargout==2:
            fwhm=np.mean(fwhm)
            return image1,fwhm 
        if nargout==1:
            return image1
         
              
    num_cores = multiprocessing.cpu_count()
    print('Number of cores = '+str(num_cores))           
    nrun=len(imgs_abs_dir[0])
    nsub=len(imgs_abs_dir)
    mask1=np.reshape(mask,[mask.shape[0]* mask.shape[1]* mask.shape[2]])
    
    list_of_arrays=[np.array(a) for a in range (0,nsub)]
    images = deepcopy(list_of_arrays)
    fwhms=np.zeros((nsub,))    
    
    results=Parallel(n_jobs=min(num_cores-2,ncore))(delayed(readImage)(img_name,mask,nrun,nargout) for img_name in imgs_abs_dir)
    
    if nargout==2:
        for i in range(0,nsub):
            images[i]=results[i][0]
            fwhms[i]=results[i][1]
            
        return images,fwhms
    
    if nargout==1:
        for i in range(0,nsub):
            images[i]=results[i]       
        
        return images
    
    


def BWAS_glm_online(image_names,targets,covariates,mask,options):
#        BWAS_glm_online(image_names,targets,covariates,mask,options)
#        Inputs: images_names: N*P cells, N subjects, P fMRI runs to connect together
#                              (within run mean and variance normalization);
#                              support nii.gz or cifti format.
#                Targets: N*Q matrix, N subjects, Q target variables,
#                         results will be stored in separate directories.
#                covariates: The covariates to regress out
#                mask: the image masks
#                      (for nii.gz, it is a 3D one, and for cifti, it is a 1D vector )
#                A list of options from previous scripts


    #where to save your results
    result_dir = options['result_dir']
    print('Results will be saved to '+result_dir)
    #subject dimension batch size
    sub_batchsize=options['sub_batchsize']
    #voxel dimension batch size
    vol_batchsize=options['vol_batchsize']
    #number of voxels
    nvol=options['nvol']
    print('Number of voxels = '+str(nvol))

    #number of runs (e.g. HCP)
    #nrun=size(image_names,2);
    #number of behavior of interests
    nsub,ntargets=targets.shape
    print('Number of Subjects = '+str(nsub))
    print('Number of Variable of Interests = '+str(ntargets))
    
    #mkdir of results
    for i in range(0,targets.shape[1]):
        if not os.path.exists(result_dir+'glm_results'+('%03d' % i)):
            os.mkdir(result_dir+'glm_results'+('%03d' % i))
        
    #get design matrix demensions
    design=np.hstack((np.reshape(targets[:,0],(nsub,1)),covariates,np.ones((nsub,1))));
    ncov=design.shape[1]
    print('Number of Covariates = '+str(ncov))

    #total number of subject batchs
    ind_end1=int(nsub/sub_batchsize-1e-4)+1;
    #total number of voxel batchs
    ind_end2=int(nvol/vol_batchsize-1e-4)+1;




    ####estimating Betas======================================================
    
    
    fname=os.path.join(result_dir,'glm_results'+ ('%03d' % (ntargets-1))+'/stat_map' + ('%04d' % (ind_end2-1)) + '_iter' + ('%03d' % (ind_end1-1)) + '_Betas.npy')    
    if os.path.exists(fname):
        print('We have got the results of all Betas ...')
        print('Omit the steps of estimating Betas...')
    else:
        
        FWHMs=np.zeros((nsub,1))        
        for i in range(0,ind_end1):
    
            st11=int(i*sub_batchsize)
            en11=min(int((i+1)*sub_batchsize),nsub)
            #sub_ids=subj_list(st11:en11);
            #load images of this subject batch
            print('Loading rsfMRI data of subject batch... '+
                    str(i+1)+'/'+str(ind_end1)+'...')
            print('Total number of subjects in this batch = '+ str(en11-st11)+'...')
    
    
            #read data, estimate smoothness
            images,fwhms=BWAS_prepare(image_names[st11:en11],mask,2)
            #print(len(images))
            FWHMs[st11:en11,0]=fwhms
            #get network of each voxel batch
            for j in range(0,ind_end2):
                start = time.time()
                #get voxel index to analysis
                st1=int(j*vol_batchsize)
                en1=min(int(j+1)*vol_batchsize,nvol)
                nn=int(nvol*(en1-st1))
                #correlation matrix
                print('Calculating the correlation matrix across all subjects... '+
                    str(j+1)+'/'+str(ind_end2)+'...')
                r1=np.zeros((en11-st11,nvol,en1-st1),dtype='float32')
                for jjj in range(0,en11-st11):
                    img_tmp=images[jjj]
                    tmp=BWAS_correlation(img_tmp,img_tmp[:,st1:en1])
                    tmp[tmp>0.9999]=0
                    r1[jjj,:,:]=tmp      
                r1=BWAS_fisher_z(r1)
                r1=np.reshape(r1,(en11-st11,nn))
                r1[np.isnan(r1)]=0
                nfc_tmp=r1.shape[1]
    
                print('Calculating the Statistical Parametirc Maps... '+
                    str(j+1)+'/'+str(ind_end2)+'...')
    
                #GLM across different behavior variables
                for ii in range(0,ntargets):
                    #get design matrix
                    design1=np.float32(np.hstack((np.reshape(targets[:,ii],(nsub,1)),covariates)))
                    
                    
                    if i==0:
    
                        #if this is not the last iter
                        if (ind_end1-1)>0:
                            
                            #This is the first batch, initialize Betas0 
                            Betas0=np.zeros((ncov,nfc_tmp),dtype='float32')
                            #Do GLM update 
                            Betas_glm=BWAS_regression_online1(design1,r1,st11,en11,Betas0)
                            fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % i) + '_Betas.npy')
                            np.save(fname1,Betas_glm)                                             
    
                    else:
                        #load results of lastest iteration
                        fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % (i-1)) + '_Betas.npy')
                        Betas0=np.load(fname1)
                        #update beta
                        Betas_glm=BWAS_regression_online1(design1,r1,st11,en11,Betas0)
                        #remove previous beta
                        os.remove(fname1)                    
                        #save results
                        fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % i) + '_Betas.npy')
                        np.save(fname1,Betas_glm)
      
                            
    
                end = time.time()
                print('Total run time for Beta = ' + str(round(end - start,2))+ ' seconds...')
            np.save(os.path.join(result_dir,'FWHMs.npy'),FWHMs)
                            

    ####estimating Residual variances======================================================

    fname=os.path.join(result_dir,'glm_results'+ ('%03d' % (ntargets-1))+'/stat_map' + ('%04d' % (ind_end2-1)) + '_iter' + ('%03d' % (ind_end1-1)) + '_Sigma.npy')    
    if os.path.exists(fname):
        print('We have got the results of all Residual Variances ...')
        print('Omit the steps of estimating Residual Variances...')
    else:
    
        for i in range(0,ind_end1):
    
            st11=int(i*sub_batchsize)
            en11=min(int((i+1)*sub_batchsize),nsub)
            #sub_ids=subj_list(st11:en11);
            #load images of this subject batch
            print('Loading rsfMRI data of subject batch... '+
                    str(i+1)+'/'+str(ind_end1)+'...')
            print('Total number of subjects in this batch = '+ str(en11-st11)+'...')
    
    
            #read data, estimate smoothness
            images=BWAS_prepare(image_names[st11:en11],mask,1)
            
            #get network of each voxel batch
            for j in range(0,ind_end2):
                start = time.time()
                #get voxel index to analysis
                st1=int(j*vol_batchsize)
                en1=min(int(j+1)*vol_batchsize,nvol)
                nn=int(nvol*(en1-st1))
                #correlation matrix
                print('Calculating the correlation matrix across all subjects... '+
                    str(j+1)+'/'+str(ind_end2)+'...')
                r1=np.zeros((en11-st11,nvol,en1-st1),dtype='float32')
                for jjj in range(0,en11-st11):
                    img_tmp=images[jjj]
                    tmp=BWAS_correlation(img_tmp,img_tmp[:,st1:en1])
                    tmp[tmp>0.9999]=0
                    r1[jjj,:,:]=tmp      
                r1=BWAS_fisher_z(r1)
                r1=np.reshape(r1,(en11-st11,nn))
                r1[np.isnan(r1)]=0
                nfc_tmp=r1.shape[1]
    
                print('Calculating the Statistical Parametirc Maps... '+
                    str(j+1)+'/'+str(ind_end2)+'...')
    
                #GLM across different behavior variables
                for ii in range(0,ntargets):
                    #get design matrix
                    design1=np.float32(np.hstack((np.reshape(targets[:,ii],(nsub,1)),covariates)))
                    
                    fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % (ind_end1-1)) + '_Betas.npy')
                    Betas_glm=np.load(fname1)
                    
                    if i==0:
    
                        #if this is not the last iter
                        if (ind_end1-1)>0:
                            
                            #This is the first batch, initialize Sigma0
                            Sigma0=np.zeros((1,nfc_tmp),dtype='float32')
                            #Do GLM update 
                            Sigma_glm=BWAS_regression_online2(design1,r1,st11,en11,Betas_glm,Sigma0)
                            fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % i) + '_Sigma.npy')
                            np.save(fname1,Sigma_glm)                                             
    
                    else:
                        #load results of lastest iteration
                        fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % (i-1)) + '_Sigma.npy')
                        Sigma0=np.load(fname1)
                        #update beta
                        Sigma_glm=BWAS_regression_online2(design1,r1,st11,en11,Betas_glm,Sigma0)
                        #remove previous beta
                        os.remove(fname1)                    
                        #save results
                        fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % i) + '_Sigma.npy')
                        np.save(fname1,Sigma_glm)
      
                            
    
                end = time.time()
                print('Total run time for Residual variance = ' + str(round(end - start,2))+ ' seconds...')

    #estimate z-stat
    for j in range(0,ind_end2):
        start = time.time()
        #get voxel index to analysis
        st1=int(j*vol_batchsize)
        en1=min(int(j+1)*vol_batchsize,nvol)

        print('Calculating the Statistical Parametirc Maps... '+
            str(j+1)+'/'+str(ind_end2)+'...')            
        #GLM across different behavior variables
        for ii in range(0,ntargets):
            #get design matrix
            design1=np.float32(np.hstack((np.reshape(targets[:,ii],(nsub,1)),covariates)))
            
            fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % (ind_end1-1)) + '_Betas.npy')
            Betas_glm=np.load(fname1)            
            fname2=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % (ind_end1-1)) + '_Sigma.npy')
            Sigma_glm=np.load(fname2)  
            os.remove(fname1)
            os.remove(fname2)
            
            stat_map=BWAS_regression_online3(design1,Betas_glm,Sigma_glm)
            #do reshape
            stat_map=np.reshape(stat_map,(nvol,(en1-st1)))
            #convert t to z
            stat_map=norm.ppf(scipy.stats.t.cdf(stat_map,nsub-ncov-1))
            #save results
            fname=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % i) + '.npy')
            np.save(fname,stat_map)            
            
        end = time.time()
        print('Total run time for z-stat = ' + str(round(end - start,2))+ ' seconds...')
    return 


##This is for parallel in CPU################################################################################################
################################################################################################################################################
def BWAS_net_parallel(img,st1,en1,nn):
    
    r=BWAS_correlation(img,img[:,st1:en1])
    r[r>0.9999]=0
    r=BWAS_fisher_z(r)
    r[np.isnan(r)]=0
    r=np.reshape(r,(nn,))
    
    return r

def BWAS_parallel_beta(j,i,result_dir,images,targets,covariates,vol_batchsize,st11,en11,nvol,nsub,ncov,ntargets,ind_end1,ind_end2):
    
    start = time.time()
    #get voxel index to analysis
    st1=int(j*vol_batchsize)
    en1=min(int(j+1)*vol_batchsize,nvol)
    nn=int(nvol*(en1-st1))
    #correlation matrix
    print('Calculating the correlation matrix across all subjects... '+
        str(j+1)+'/'+str(ind_end2)+'...')
    r1=np.zeros((en11-st11,nvol,en1-st1),dtype='float32')
    for jjj in range(0,en11-st11):
        img_tmp=images[jjj]
        tmp=BWAS_correlation(img_tmp,img_tmp[:,st1:en1])
        tmp[tmp>0.9999]=0
        r1[jjj,:,:]=tmp      
    r1=BWAS_fisher_z(r1)
    r1=np.reshape(r1,(en11-st11,nn))
    r1[np.isnan(r1)]=0
    nfc_tmp=r1.shape[1]

    print('Calculating the Statistical Parametirc Maps... '+
        str(j+1)+'/'+str(ind_end2)+'...')

    #GLM across different behavior variables
    for ii in range(0,ntargets):
        #get design matrix
        design1=np.float32(np.hstack((np.reshape(targets[:,ii],(nsub,1)),covariates)))
        
        
        if i==0:

            #if this is not the last iter
            if (ind_end1-1)>0:
                
                #This is the first batch, initialize Betas0 
                Betas0=np.zeros((ncov,nfc_tmp),dtype='float32')
                #Do GLM update 
                Betas_glm=BWAS_regression_online1(design1,r1,st11,en11,Betas0)
                fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % i) + '_Betas.npy')
                np.save(fname1,Betas_glm)                                             

        else:
            #load results of lastest iteration
            fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % (i-1)) + '_Betas.npy')
            Betas0=np.load(fname1)
            #update beta
            Betas_glm=BWAS_regression_online1(design1,r1,st11,en11,Betas0)
            #remove previous beta
            os.remove(fname1)                    
            #save results
            fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % i) + '_Betas.npy')
            np.save(fname1,Betas_glm)
  
                

    end = time.time()
    print('Total run time for Beta = ' + str(round(end - start,2))+ ' seconds...')


def BWAS_parallel_Sigma(j,i,result_dir,images,targets,covariates,vol_batchsize,st11,en11,nvol,nsub,ncov,ntargets,ind_end1,ind_end2):
    
    
    start = time.time()
    #get voxel index to analysis
    st1=int(j*vol_batchsize)
    en1=min(int(j+1)*vol_batchsize,nvol)
    nn=int(nvol*(en1-st1))
    #correlation matrix
    print('Calculating the correlation matrix across all subjects... '+
        str(j+1)+'/'+str(ind_end2)+'...')
    r1=np.zeros((en11-st11,nvol,en1-st1),dtype='float32')
    for jjj in range(0,en11-st11):
        img_tmp=images[jjj]
        tmp=BWAS_correlation(img_tmp,img_tmp[:,st1:en1])
        tmp[tmp>0.9999]=0
        r1[jjj,:,:]=tmp      
    r1=BWAS_fisher_z(r1)
    r1=np.reshape(r1,(en11-st11,nn))
    r1[np.isnan(r1)]=0
    nfc_tmp=r1.shape[1]
    

    print('Calculating the Statistical Parametirc Maps... '+
        str(j+1)+'/'+str(ind_end2)+'...')

    #GLM across different behavior variables
    for ii in range(0,ntargets):
        #get design matrix
        design1=np.float32(np.hstack((np.reshape(targets[:,ii],(nsub,1)),covariates)))
        
        fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % (ind_end1-1)) + '_Betas.npy')
        Betas_glm=np.load(fname1)
        
        if i==0:

            #if this is not the last iter
            if (ind_end1-1)>0:
                
                #This is the first batch, initialize Sigma0
                Sigma0=np.zeros((1,nfc_tmp),dtype='float32')
                #Do GLM update 
                Sigma_glm=BWAS_regression_online2(design1,r1,st11,en11,Betas_glm,Sigma0)
                fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % i) + '_Sigma.npy')
                np.save(fname1,Sigma_glm)                                             

        else:
            #load results of lastest iteration
            fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % (i-1)) + '_Sigma.npy')
            Sigma0=np.load(fname1)
            #update beta
            Sigma_glm=BWAS_regression_online2(design1,r1,st11,en11,Betas_glm,Sigma0)
            #remove previous beta
            os.remove(fname1)                    
            #save results
            fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % i) + '_Sigma.npy')
            np.save(fname1,Sigma_glm)
                  
    end = time.time()
    print('Total run time for Residual variance = ' + str(round(end - start,2))+ ' seconds...')



def BWAS_parallel_zstat(j,i,result_dir,targets,covariates,vol_batchsize,st11,en11,nvol,nsub,ncov,ntargets,ind_end1,ind_end2):
    
    start = time.time()
    #get voxel index to analysis
    st1=int(j*vol_batchsize)
    en1=min(int(j+1)*vol_batchsize,nvol)

    print('Calculating the Statistical Parametirc Maps... '+
        str(j+1)+'/'+str(ind_end2)+'...')            
    #GLM across different behavior variables
    for ii in range(0,ntargets):
        #get design matrix
        design1=np.float32(np.hstack((np.reshape(targets[:,ii],(nsub,1)),covariates)))
        
        fname1=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % (ind_end1-1)) + '_Betas.npy')
        Betas_glm=np.load(fname1)            
        fname2=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % (ind_end1-1)) + '_Sigma.npy')
        Sigma_glm=np.load(fname2)  
        os.remove(fname1)
        os.remove(fname2)
        
        stat_map=BWAS_regression_online3(design1,Betas_glm,Sigma_glm)
        #do reshape
        stat_map=np.reshape(stat_map,(nvol,(en1-st1)))
        #convert t to z
        stat_map=norm.ppf(scipy.stats.t.cdf(stat_map,nsub-ncov-1))
        #save results
        fname=os.path.join(result_dir,'glm_results'+ ('%03d' % ii)+'/stat_map' + ('%04d' % j) + '_iter' + ('%03d' % i) + '.npy')
        np.save(fname,stat_map)            
        
    end = time.time()
    print('Total run time for z-stat = ' + str(round(end - start,2))+ ' seconds...')

    return


def BWAS_glm_online_parallel(image_names,targets,covariates,mask,options):
#        BWAS_glm_online(image_names,targets,covariates,mask,options)
#        Inputs: images_names: N*P cells, N subjects, P fMRI runs to connect together
#                              (within run mean and variance normalization);
#                              support nii.gz or cifti format.
#                Targets: N*Q matrix, N subjects, Q target variables,
#                         results will be stored in separate directories.
#                covariates: The covariates to regress out
#                mask: the image masks
#                      (for nii.gz, it is a 3D one, and for cifti, it is a 1D vector )
#                A list of options from previous scripts


    #where to save your results
    result_dir = options['result_dir']
    print('Results will be saved to '+result_dir)
    #subject dimension batch size
    sub_batchsize=options['sub_batchsize']
    #voxel dimension batch size
    vol_batchsize=options['vol_batchsize']
    #number of voxels
    nvol=options['nvol']
    print('Number of voxels = '+str(nvol))

    #number of runs (e.g. HCP)
    #nrun=size(image_names,2);
    #number of behavior of interests
    nsub,ntargets=targets.shape
    print('Number of Subjects = '+str(nsub))
    print('Number of Variable of Interests = '+str(ntargets))
    
    #mkdir of results
    for i in range(0,targets.shape[1]):
        if not os.path.exists(result_dir+'glm_results'+('%03d' % i)):
            os.mkdir(result_dir+'glm_results'+('%03d' % i))
        
    #get design matrix demensions
    design=np.hstack((np.reshape(targets[:,0],(nsub,1)),covariates));
    ncov=design.shape[1]
    print('Number of Covariates (including a column of ones) = '+str(ncov))

    #total number of subject batchs
    ind_end1=int(nsub/sub_batchsize-1e-4)+1;
    #total number of voxel batchs
    ind_end2=int(nvol/vol_batchsize-1e-4)+1;




    ####estimating Betas======================================================
    
    
    fname=os.path.join(result_dir,'glm_results'+ ('%03d' % (ntargets-1))+'/stat_map' + ('%04d' % (ind_end2-1)) + '_iter' + ('%03d' % (ind_end1-1)) + '_Betas.npy')    
    if os.path.exists(fname):
        print('We have got the results of all Betas ...')
        print('Omit the steps of estimating Betas...')
    else:
        
        FWHMs=np.zeros((nsub,1))        
        for i in range(0,ind_end1):
    
            st11=int(i*sub_batchsize)
            en11=min(int((i+1)*sub_batchsize),nsub)
            #sub_ids=subj_list(st11:en11);
            #load images of this subject batch
            print('Loading rsfMRI data of subject batch... '+
                    str(i+1)+'/'+str(ind_end1)+'...')
            print('Total number of subjects in this batch = '+ str(en11-st11)+'...')
    
    
            #read data, estimate smoothness
            images,fwhms=BWAS_prepare_parallel(image_names[st11:en11],mask,options['ncore'],2)
            #print(len(images))
            FWHMs[st11:en11,0]=fwhms
            #get network of each voxel batch
            #for j in range(0,ind_end2):
            num_cores = multiprocessing.cpu_count()
            Parallel(n_jobs=min(num_cores-2,options['ncore']))(delayed(
                    BWAS_parallel_beta)(j,i,result_dir,images,targets,covariates,vol_batchsize,st11,en11,nvol,nsub,ncov,ntargets,ind_end1,ind_end2) for j in range(0,ind_end2))
            
            
            np.save(os.path.join(result_dir,'FWHMs.npy'),FWHMs)
                            

    ####estimating Residual variances======================================================

    fname=os.path.join(result_dir,'glm_results'+ ('%03d' % (ntargets-1))+'/stat_map' + ('%04d' % (ind_end2-1)) + '_iter' + ('%03d' % (ind_end1-1)) + '_Sigma.npy')    
    if os.path.exists(fname):
        print('We have got the results of all Residual Variances ...')
        print('Omit the steps of estimating Residual Variances...')
    else:
    
        for i in range(0,ind_end1):
    
            st11=int(i*sub_batchsize)
            en11=min(int((i+1)*sub_batchsize),nsub)
            #sub_ids=subj_list(st11:en11);
            #load images of this subject batch
            print('Loading rsfMRI data of subject batch... '+
                    str(i+1)+'/'+str(ind_end1)+'...')
            print('Total number of subjects in this batch = '+ str(en11-st11)+'...')
                
            #read data, estimate smoothness
            images=BWAS_prepare_parallel(image_names[st11:en11],mask,options['ncore'],1)
            #get network of each voxel batch
            #for j in range(0,ind_end2):
            num_cores = multiprocessing.cpu_count()
            Parallel(n_jobs=min(num_cores-2,options['ncore']))(delayed(
                    BWAS_parallel_Sigma)(j,i,result_dir,images,targets,covariates,vol_batchsize,st11,en11,nvol,nsub,ncov,ntargets,ind_end1,ind_end2) for j in range(0,ind_end2))
            
            
            
    #estimate z-stat
    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=min(num_cores-2,options['ncore']))(delayed(
            BWAS_parallel_zstat)(j,i,result_dir,targets,covariates,vol_batchsize,st11,en11,nvol,nsub,ncov,ntargets,ind_end1,ind_end2) for j in range(0,ind_end2))
    

##This is the end for parallel in CPU  1################################################


    
    

def BWAS_GRF_6D_density(t):
    #This is EC density of 6D Gaussian Random Field
    
    EC=np.zeros((7,1))
    
    a = 4.0 * np.log(2.0)
    a1= 2.0 * np.pi
    b = np.exp(-np.square(t)/2.0)
    
    EC[0] = 1.0- norm.cdf(t) 
    EC[1] = a ** (1.0/2.0) / a1 * b 
    EC[2] = a/(a1 ** (3.0/2.0)) * b * t 
    EC[3] = a ** (3.0/2.0)/(a1**2.0) * b * (np.power(t,2) - 1)
    EC[4] = a ** 2.0/(a1 ** (5.0/2.0)) * b * ( np.power(t,3) - 3 * t )
    EC[5] = a ** (5.0/2.0)/(a1 ** 3.0) * b * ( np.power(t,4) - 6 * np.power(t,2) + 3 )
    EC[6] = a ** 3.0/(a1 ** (7.0/2.0)) * b * ( np.power(t,5) - 10 * np.power(t,3) + 15*t)
    
    return EC


def BWAS_whole_brain_cluster_p(nvoxel1,nvoxel2,cdt_z,cl_size,fwhm):
    
    fwhm=np.array(fwhm,dtype='float64')
    cdt_z=np.array(cdt_z,dtype='float64')
    cl_size=np.array(cl_size,dtype='float64')
    
    r1=(nvoxel1/(4.0/3.0*np.pi))**(1/3.0)    
    R1=[1, 
        4*(r1**3/np.prod(fwhm))**(1/3.0),
        2*np.pi*(r1**3/np.prod(fwhm))**(2/3.0), 
        nvoxel1/np.prod(fwhm)]
        
    r2=(nvoxel2/(4.0/3.0*np.pi))**(1/3.0)    
    R2=[1, 
        4*(r2**3/np.prod(fwhm))**(1/3.0),
        2*np.pi*(r2**3/np.prod(fwhm))**(2/3.0), 
        nvoxel2/np.prod(fwhm)]    
        
    V=nvoxel1*nvoxel2
    
    EC=BWAS_GRF_6D_density(cdt_z);
     
     
    EN=V*(norm.cdf(-np.abs(cdt_z)))
    
    EL1=0
    for i in range(0,4):
        for j in range(0,4):
            EL1=EL1+R1[i]*R2[j]*EC[i+j];

    
    if EL1<0:
        p_uncorrected=1
        p_corrected=1
    else:
        phi=(scipy.special.gamma(6.0/2.0+1)*EL1/EN)**(2.0/6.0);
        
        p_uncorrected=np.exp(-phi*cl_size**(2.0/6.0));
        
        p_corrected=1-np.exp(-EL1*p_uncorrected)
 
    return p_corrected,p_uncorrected

  
def BWAS_get_sparse_dismat(dims):
    #dims n*p
    #output n*n
    
    nvol=dims.shape[0]
    nstep=10000
    
    for i in range(0,nvol,nstep):
        s=i
        s1=min(i+nstep,nvol)
        
        tmptmp=cdist(dims[s:s1,:],dims,'euclidean')
        #tmp=np.where((tmptmp<1.5) & (tmptmp>0))
        tmp=np.where(tmptmp<1.5)
        
        if i==0:
            coo1=tmp[0]
            coo2=tmp[1]
        else:            
            coo1=np.append(coo1,tmp[0]+i)
            coo2=np.append(coo2,tmp[1])
        #print(i)

    adj=sparse.csc_matrix((np.ones((len(coo1))), (coo1, coo2)), shape=(nvol, nvol))    
    
    return adj
    
   
    


def BWAS_analysis_result(result_dir,mask_file,CDT,fwhm,options):
    
    #result_dir='/Users/wgong/Documents/MATLAB/hcp_bwas/glm_results004/'
    #iter_number=2
    #mask_file='/Users/wgong/Documents/MATLAB/hcp_bwas/aal2_4mm.nii.gz'
    #CDT=4.5
    #fwhm=[2,2,2]
    
    info=nib.load(mask_file)
    mask=np.float64(info.get_data()!=0)
    affine=info.affine

    #mask=np.float32((mask>0) & (mask<9000))
    #img = nib.Nifti1Image(mask, affine)
    #nib.save(img, result_dir+'hcp_bwas_mask_final_4mm.nii.gz')
      
    dim=np.array(np.where(mask!=0)).T
    print('Number of voxels = '+str(dim.shape[0]))

    mni=np.dot(np.hstack((dim,np.ones((dim.shape[0],1)))),affine.T[:,0:3])
    
    nfile=len(glob.glob(result_dir+'stat_map*.npy'))
    files_to_read=sorted(glob.glob(result_dir+'stat_map*.npy'))
    print('Number of parts = '+str(nfile))
    
    #load all fc z-value
    print('Reading stat maps...')
    ind_sig=np.zeros((dim.shape[0],dim.shape[0]),dtype='uint8')
    
    nvol=dim.shape[0]
    vol_batchsize=options['vol_batchsize']
    
    for i in range(0,nfile):
        st1=int(i*vol_batchsize)
        en1=min(int((i+1)*vol_batchsize),int(nvol))
        
        stat_map=np.load(files_to_read[i])
        stat_map[np.isnan(stat_map)]=0
        #significant
        tmp=np.abs(stat_map)>CDT
        ind_sig[:,st1:en1]=tmp
        
        coo_tmp=np.array(np.where(tmp==1)).T
        coo_tmp[:,1]=coo_tmp[:,1]+st1
        zval_tmp=np.array(stat_map[tmp],ndmin=2).T
        data_tmp=np.hstack((coo_tmp,zval_tmp))
        
        if i==0:
            data=data_tmp            
        else:
            data=np.vstack((data,data_tmp))
            
        print('Loading and preprocessing '+files_to_read[i])
    
    print('Done...')
    stat_map=sparse.csc_matrix((data[:,2], (data[:,0], data[:,1])), shape=(nvol, nvol))
    
 
    #index of fcs with |z|>cdt
    
    ind_sig=np.triu(ind_sig,0)
    #matrix index
    ind_sig_mat=np.array(np.where(ind_sig==1)).T
    
    ind_sig=None
        
    
    print('Number of FCs with |z|>'+str(CDT)+' = '+str(ind_sig_mat.shape[0]))
    
    if ind_sig_mat.shape[0]==0:
        print('No significant FCs...')
    else:
        
        #get 6d coordiantes
        FC_6d_coordiantes=np.hstack((dim[ind_sig_mat[:,0],:],dim[ind_sig_mat[:,1],:]))
        FC_6d_MNI=np.hstack((mni[ind_sig_mat[:,0],:],mni[ind_sig_mat[:,1],:]))

    
        #calculate (sparse) distance matrix of mask
        adj=BWAS_get_sparse_dismat(dim)
        
        #coo=np.where(cdist(dim,dim,'euclidean')<1.5)
        #coo=np.array(coo).T
        #adj=sparse.csc_matrix((np.ones((len(coo))), (coo[:,0], coo[:,1])), shape=(nvol, nvol))
        
        #first 3 dims' adj matrix
        adj1=adj[ind_sig_mat[:,0]].T[ind_sig_mat[:,0]]
        adj2=adj[ind_sig_mat[:,1]].T[ind_sig_mat[:,1]]
                
        #get final intersected adj mat
        adj_final=(adj1+adj2)==2
        
        #get connected components
        n_components,labels=sparse.csgraph.connected_components(adj_final,directed=False)
        print('Number of FC clusters = '+str(n_components))
        
        
        result_table=np.zeros((n_components,6))
        table_title=['Number of FC','pval corrected','pval uncorrected',
                     'max zvalue','Number of voxels in region1',
                     'Number of voxels in region2']
        infos_tmp=[]
        tmptmp=[];
        nn=0;
        for j in range(0,n_components):
            #print(j)
            nFC=np.sum(labels==j)
            p_corrected,p_uncorrected=BWAS_whole_brain_cluster_p(dim.shape[0]/np.sqrt(2),dim.shape[0]/np.sqrt(2),CDT,nFC,fwhm)
            
            index1=ind_sig_mat[labels==j,0]
            index2=ind_sig_mat[labels==j,1]
            FC_6d_coordiantes_cl=FC_6d_coordiantes[labels==j,:]
            zval=np.matrix(stat_map[index1,index2]).T
            
            
            
            if p_corrected<0.05:
                infos_tmp.append(np.hstack((FC_6d_coordiantes_cl,zval)))
                tmptmp.append(np.array(nFC))
                
                nn=nn+1
                if nn==1:
                    FC_6d_MNI_cl=np.hstack((FC_6d_MNI[labels==j,0:6],zval))
                    FC_6d_mat_cl=np.hstack((FC_6d_coordiantes[labels==j,0:6],zval))
                else:
                    FC_6d_MNI_cl=np.vstack((FC_6d_MNI_cl,np.hstack((FC_6d_MNI[labels==j,0:6],zval))))
                    FC_6d_mat_cl=np.vstack((FC_6d_mat_cl,np.hstack((FC_6d_coordiantes[labels==j,0:6],zval))))

            nvol1=len(np.unique(index1))
            nvol2=len(np.unique(index2))
            result_table[j,0]=nFC
            result_table[j,1]=p_corrected
            result_table[j,2]=p_uncorrected
            result_table[j,3]=np.max(zval)
            result_table[j,4]=nvol1
            result_table[j,5]=nvol2
            
        print('Number of Significant FC clusters = '+str(len(infos_tmp)))

            
        ind=np.argsort(result_table[:,2])
        ind1=np.argsort(-np.array(tmptmp))
        
        result_table=result_table[ind,:]
        
        infos=[]        
        for j in range(0,len(ind1)):
            infos.append(infos_tmp[ind1[j]])
            print(infos_tmp[ind1[j]].shape)
        
        ma_3d_cl=np.zeros((mask.shape[0],mask.shape[1],mask.shape[2]))
        
        if np.sum(result_table[:,1]<0.05)>0:
            FC_6d_mat_cl=np.int64(FC_6d_mat_cl)
            for j in range(0,FC_6d_mat_cl.shape[0]):
                ma_3d_cl[FC_6d_mat_cl[j,0],FC_6d_mat_cl[j,1],FC_6d_mat_cl[j,2]]=ma_3d_cl[FC_6d_mat_cl[j,0],FC_6d_mat_cl[j,1],FC_6d_mat_cl[j,2]]+1
                ma_3d_cl[FC_6d_mat_cl[j,3],FC_6d_mat_cl[j,4],FC_6d_mat_cl[j,5]]=ma_3d_cl[FC_6d_mat_cl[j,3],FC_6d_mat_cl[j,4],FC_6d_mat_cl[j,5]]+1
        else:
            FC_6d_MNI_cl=0
            FC_6d_mat_cl=0
           
        img = nib.Nifti1Image(ma_3d_cl, affine)
        nib.save(img, result_dir+'MA_CDT='+str(CDT)+'.nii.gz')
  
        
        np.savez(result_dir+'BWAS_result_CDT='+str(CDT)+'.npz',infos=infos,result_table=result_table,table_title=table_title,FC_6d_MNI_cl=FC_6d_MNI_cl,FC_6d_mat_cl=FC_6d_mat_cl)
        sio.savemat(result_dir+'BWAS_result_CDT='+str(CDT)+'.mat',{'infos':infos,'result_table':result_table,'table_title':table_title,'FC_6d_MNI_cl':FC_6d_MNI_cl,'FC_6d_mat_cl':FC_6d_mat_cl})
        
        #plot MA map
        plotting.plot_stat_map(result_dir+'MA_CDT='+str(CDT)+'.nii.gz',threshold=0,display_mode="z", 
                       cut_coords=10,colorbar=True,
                       title='Number of Significant FCs per voxel',
                       output_file=result_dir+'MA_CDT='+str(CDT)+'_z.jpg')
        plotting.plot_stat_map(result_dir+'MA_CDT='+str(CDT)+'.nii.gz',threshold=0,display_mode="x", 
                       cut_coords=10,colorbar=True,
                       title='Number of Significant FCs per voxel',
                       output_file=result_dir+'MA_CDT='+str(CDT)+'_x.jpg')
        plotting.plot_stat_map(result_dir+'MA_CDT='+str(CDT)+'.nii.gz',threshold=0,display_mode="y", 
                       cut_coords=10,colorbar=True,
                       title='Number of Significant FCs per voxel',
                       output_file=result_dir+'MA_CDT='+str(CDT)+'_y.jpg')
        
        toolbox_path=options['toolbox_path']
        #BWAS_plot_FC_clusters
        if np.sum(result_table[:,1]<0.05)>0:
            if not os.path.exists(result_dir+'FC_figure'):
                os.mkdir(result_dir+'FC_figure')
                os.mkdir(result_dir+'FC_cluster_files')
                
            for i in range(0,len(infos)):
            
                stat_map1=np.zeros((mask.shape[0],mask.shape[1],mask.shape[2]))
                for j in range(0,infos[i].shape[0]):
                    stat_map1[int(infos[i][j,0]),int(infos[i][j,1]),int(infos[i][j,2])]=stat_map1[int(infos[i][j,0]),int(infos[i][j,1]),int(infos[i][j,2])]+1*(np.sign(infos[i][j,6]))
                img = nib.Nifti1Image(stat_map1, affine)
                nib.save(img, result_dir+'img_to_plot.nii.gz')                
                fig=plotting.plot_stat_map(result_dir+'img_to_plot.nii.gz',bg_img=toolbox_path+'MNI152_T1_1mm_brain.nii.gz',threshold=0,display_mode="z", 
                       cut_coords=1,colorbar=False,black_bg=False)                  
                fig.title('FC Cluster '+str(i+1)+', nFC = '+str(infos[i].shape[0])+', FWER p-value = '+str(np.format_float_scientific(result_table[i,1],2)),size=5)
                fig.savefig(result_dir+'1.jpg',dpi=512)  
                     
                img1 = np.array(misc.imread(result_dir+'1.jpg'),dtype='uint8')
                os.rename(result_dir+'img_to_plot.nii.gz',result_dir+'/FC_cluster_files/FC_cluster_'+('%04d' % (i+1))+'_1.nii.gz')
                
                stat_map1=np.zeros((mask.shape[0],mask.shape[1],mask.shape[2]))
                for j in range(0,infos[i].shape[0]):
                    stat_map1[int(infos[i][j,3]),int(infos[i][j,4]),int(infos[i][j,5])]=stat_map1[int(infos[i][j,3]),int(infos[i][j,4]),int(infos[i][j,5])]+1*(np.sign(infos[i][j,6]))
                    
                img = nib.Nifti1Image(stat_map1, affine)
                nib.save(img, result_dir+'img_to_plot.nii.gz')
                fig=plotting.plot_stat_map(result_dir+'img_to_plot.nii.gz',bg_img=toolbox_path+'MNI152_T1_1mm_brain.nii.gz',threshold=0,display_mode="z", 
                       cut_coords=1,colorbar=True,black_bg=False)  
                fig.savefig(result_dir+'1.jpg',dpi=512)  
                
                img2 = np.array(misc.imread(result_dir+'1.jpg'),dtype='uint8')
                os.rename(result_dir+'img_to_plot.nii.gz',result_dir+'/FC_cluster_files/FC_cluster_'+('%04d' % (i+1))+'_2.nii.gz')
                
                
                if i==0:
                    img_tmp=np.hstack((img1,img2))
                    #image_final=img_tmp
                    misc.imsave(result_dir+'FC_figure/BWAS_FC_clusters_'+('%04d' % (i+1))+'.pdf',img_tmp) 
                else:
                    img_tmp=np.hstack((img1,img2))
                    #image_final=np.vstack((image_final,img_tmp))
                    misc.imsave(result_dir+'FC_figure/BWAS_FC_clusters_'+('%04d' % (i+1))+'.pdf',img_tmp) 
                    
            #misc.imsave(result_dir+'BWAS_FC_clusters_plot.jpg',image_final)     
            merger = PdfFileMerger()
            for i in range(0,len(infos)):
                merger.append(open(result_dir+'FC_figure/BWAS_FC_clusters_'+('%04d' % (i+1))+'.pdf', 'rb'))

            with open(result_dir+'BWAS_FC_clusters_plot.pdf', 'wb') as fout:
                merger.write(fout)  
            os.remove(result_dir+'1.jpg')  
            #os.remove(result_dir+'img_to_plot.nii.gz')
            
        #for braingl only
        if np.sum(result_table[:,1]<0.05)>0:
            FC_6d_MNI_cl1=FC_6d_MNI_cl+[93,125,75,93,125,75,0]
            ind11=np.where(FC_6d_MNI_cl1[:,6]<0)[0]
            np.savetxt(result_dir+'brain_gl_negative.txt',FC_6d_MNI_cl1[ind11,:],
                       delimiter=' ')
            os.rename(result_dir+'brain_gl_negative.txt',result_dir+'brain_gl_negative.cxls')
            
            FC_6d_MNI_cl1=FC_6d_MNI_cl+[93,125,75,93,125,75,0]
            ind11=np.where(FC_6d_MNI_cl1[:,6]>0)[0]
            np.savetxt(result_dir+'brain_gl_positive.txt',FC_6d_MNI_cl1[ind11,:],
                       delimiter=' ')
            os.rename(result_dir+'brain_gl_positive.txt',result_dir+'brain_gl_positive.cxls')
            
        print('All analysis finished!')
        
        
        
    return



  

def BWAS_run_full_analysis(result_dir,image_dir,mask_file,toolbox_path,targets_file,cov_file,CDT,memory_limit_per_core=16,ncore=1):
    

    result_dir=os.path.join(result_dir,'')
    
    #get image file names
    #sort them
    if isinstance(image_dir,str):
        image_dir=os.path.join(image_dir,'')
        
        img_file_names=sorted(glob.glob(image_dir+'*.nii.gz'))
        image_names=[]
        for i in range(0,len(img_file_names)):
            tmp=[]
            tmp.append(img_file_names[i])        
            image_names.append(tmp)
    else:
        image_names=image_dir   
    
    if '.npy' in targets_file:
        targets=np.load(targets_file)
        if len(targets.shape)==1:
            targets=np.array(targets,ndmin=2).T
    elif '.txt' in targets_file:
        targets=np.loadtxt(targets_file)
        if len(targets.shape)==1:
            targets=np.array(targets,ndmin=2).T
    else:
        print('Wrong "targets" file type...')
    
    if cov_file==None:
        covariates=np.ones((targets.shape[0],1))
    else:
        if '.npy' in targets_file:
            covs_tmp=np.load(cov_file)
            if len(covs_tmp.shape)==1:
                covs_tmp=np.array(covs_tmp,ndmin=2).T
            covariates=np.hstack((covs_tmp,np.ones((targets.shape[0],1))))
        elif '.txt' in targets_file:
            covs_tmp=np.loadtxt(cov_file)
            if len(covs_tmp.shape)==1:
                covs_tmp=np.array(covs_tmp,ndmin=2).T            
            covariates=np.hstack((covs_tmp,np.ones((targets.shape[0],1))))
        else:
            print('Wrong "covariates" file type...')   
            
            
    #load the mask for future use
    info=nib.load(mask_file)
    mask=info.get_data()
    mask=np.float64(mask!=0)
    nvol=int(np.sum(np.sum(mask)))
    print('Number of Voxels in mask = '+str(nvol))
    print('Shape of the Mask = '+str(mask.shape))

    #load a sample image
    info=nib.load(image_names[0][0])
    fmri=BWAS_prepare_parallel(image_names[0:1],mask,1,1)
    print(fmri[0].shape)
    gb1sub=sys.getsizeof(fmri[0])
    
    vol_batchsize=200
    sub_batchsize=min(int(memory_limit_per_core/((gb1sub+nvol*vol_batchsize*(1+covariates.shape[1])*8)/(1024.0**3))),targets.shape[0]/2)
        
    
    options={'CDT': CDT, 'sub_batchsize':sub_batchsize, 'nvol': nvol ,'vol_batchsize':vol_batchsize,
        'result_dir':result_dir,'ncore':ncore,'toolbox_path':toolbox_path}
    
    
    #if the glm analysis is not finished yet, rerun all glm,
    #else only compute the stats
    #total number of subject batchs
    ind_end1=int(targets.shape[0]/sub_batchsize-1e-4)+1;
    #total number of voxel batchs
    ind_end2=int(nvol/vol_batchsize-1e-4)+1;
    
    fname=os.path.join(result_dir,'glm_results'+ ('%03d' % (targets.shape[1]-1))+'/stat_map' + ('%04d' % (ind_end2-1)) + '_iter' + ('%03d' % (ind_end1-1)) + '.npy')
        
    if os.path.exists(fname)==0:
        if ncore==1:
            BWAS_glm_online(image_names,targets,covariates,mask,options)
        else:
            BWAS_glm_online_parallel(image_names,targets,covariates,mask,options)
        
    
    
    fwhm1=np.load(os.path.join(result_dir,'FWHMs.npy'))
    fwhm2=np.mean(fwhm1)
    if fwhm2<2.0:
        fwhm2=np.array([2.0])
    fwhm=[fwhm2,fwhm2,fwhm2]
    
    folderdir=[]
    for j in range(0,targets.shape[1]):   
        folderdir.append(result_dir+'glm_results'+('%03d' % j)+'/')
    
    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=min(num_cores-2,options['ncore']))(delayed(
            BWAS_analysis_result)(folderdir[j],mask_file,CDT,fwhm,options) for j in range(0,targets.shape[1]))
            
    #BWAS_analysis_result(folderdir,mask_file,CDT,fwhm,options)
    
    return







    
    
    
    
    
    