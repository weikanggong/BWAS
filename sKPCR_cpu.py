
# Script name: sKPCR_cpu.py
#
# Description: Functions to run sKPCR analysis using CPU
#
# Author: Weikang Gong
#
#Gong, Weikang, et al. "A powerful and efficient multivariate approach for voxel-level connectome-wide association studies." NeuroImage (2018).
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




import numpy as np

import nibabel as nib


from scipy.spatial.distance import cdist
from scipy import sparse
import scipy.sparse.linalg
from scipy.stats import rankdata
import scipy

import time

import os
import glob

from copy import deepcopy

from joblib import Parallel, delayed
import multiprocessing


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
                print(image1.shape)
                
                print('Done...')
           
            images[i]=image1
            print(images[i].shape[1])
            
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
    


def sKPCR_prepare_Laplacian(mask):
    #mask:   a 3D 0-1 mask for nii data
    #output: the laplacian matrix
    dim=np.where(mask!=0)
    dims=np.vstack((dim[0],dim[1],dim[2]))
    nvol=len(dims[0,:])
    
    nstep=10000
    for i in range(0,nvol,nstep):
        s=i
        s1=min(i+nstep,nvol)
        
        dist_map=cdist(dims[:,s:s1].T,dims.T,'euclidean');
        tmp=np.where((dist_map<1.01) & (dist_map>0))
        if i==0:
            coo1=tmp[0]
            coo2=tmp[1]
        else:            
            coo1=np.append(coo1,tmp[0]+i)
            coo2=np.append(coo2,tmp[1])
        #print(i)

    adj=sparse.csc_matrix((np.ones((len(coo1))), (coo1, coo2)), shape=(nvol, nvol))
    
    laplacian_mat=sparse.csgraph.laplacian(adj)
    
    return laplacian_mat
    

def sKPCR_linear(X,laplacian_mat,K):
    #X is n*p, laplacian_mat is p*p sparse, k is a number
    n=X.shape[0]
    Kernel0=np.dot(X*laplacian_mat,X.T)
    oneN=np.divide(np.ones((n,n)),np.float64(n))
    Kernel=np.divide((Kernel0-np.dot(oneN,Kernel0)-np.dot(Kernel0,oneN)+np.dot(np.dot(oneN,Kernel0),oneN)),np.float64(n))

    #D,U=scipy.sparse.linalg.eigs(Kernel,k=K)
    #D,U=np.linalg.eig(Kernel)
    
    D,U=scipy.linalg.eigh(Kernel,eigvals=(n-K,n-1))
    D=np.real(D)
    U=np.real(U)
    indx1=np.argsort(-D)
    D=D[indx1]
    U=U[:,indx1]
    
    
    U=np.array(np.real(U[:,0:K]),dtype='float64')
      
    return U
    

def sKPCR_regression(X,Y,cov):

    contrast=np.transpose(np.hstack(  ( np.eye(X.shape[1],X.shape[1]) , np.zeros((X.shape[1],cov.shape[1])) ))   )
    contrast=np.array(contrast,dtype='float32')
    
    design=np.hstack((X,cov))
    #degree of freedom
    df=design.shape[0]-design.shape[1]
    #
    ss=np.linalg.inv(np.dot(np.transpose(design),design))

    beta=np.dot(np.dot(ss,np.transpose(design)),Y)

    Res=Y-np.dot(design,beta)

    sigma=np.reshape(np.sqrt(np.divide(np.sum(np.square(Res),axis=0),df)),(1,beta.shape[1]))

    tmp1=np.dot(beta.T,contrast)
    tmp2=np.array(np.diag(np.dot(np.dot(contrast.T,ss),contrast)),ndmin=2)

    Tstat=np.divide(tmp1,np.dot(sigma.T,np.sqrt(tmp2)  ))


    return Tstat


def BWAS_regression(X,Y):
    #contrast of interst
    contrast=np.transpose(np.hstack((np.ones((1,1)),np.zeros((1,X.shape[1]-1)))))
    #contrast=np.array(contrast,dtype='float32')
    #degree of freedom
    df=np.float32(X.shape[0]-X.shape[1])
    #
    ss=np.linalg.pinv(np.dot(np.transpose(X),X))
    
    beta=np.dot(np.dot(ss,np.transpose(X)),Y)
    
    Res=Y-np.dot(X,beta)
      
    sigma=np.reshape(np.sqrt(np.divide(np.sum(np.square(Res),axis=0),df)),(1,beta.shape[1]))
    
    Tstat=np.divide(np.dot(np.transpose(beta),contrast),np.dot(np.transpose(sigma),np.sqrt(np.dot(np.dot(np.transpose(contrast),ss),contrast))))

    return Tstat


def sKPCR_adptive_regression(pcs,pheno,cov,numperms):
    N,P=pcs.shape
    #pheno the first one is targets
    #the 2 to end are permuted
    
    if pheno.shape[1]==1:
        yperms = np.zeros((N,numperms),dtype='float32')
        yperms[:,0]=pheno.flatten()
        for j in range(1,numperms):
            yperms[:,j] = pheno[np.random.permutation(N),0]
    else:
        yperms=pheno*1.0
        numperms=yperms.shape[1]
        
    Tstat2=np.zeros((numperms,P),dtype='float32')
    '''
    for j in range(0,P):
        design=np.hstack((np.array(pcs[:,j],ndmin=2).T,cov))
        design=np.array(design,dtype='float32')
        Tstat2[:,j]=np.square(BWAS_regression(design,yperms)).flatten()
    '''
    
    for j in range(0,numperms):
        design=np.hstack((np.array(yperms[:,j],ndmin=2).T,cov))
        design=np.array(design,dtype='float32')
        Tstat2[j,:]=np.square(BWAS_regression(design,pcs)).flatten()
    
    
    pPerm0=np.ones((P,1),dtype='float32') 
    T2stats=np.ones((P,),dtype='float32') 
    for j in range(0,P):
        #test statistic (1+numperm)*1 vector
        T0s = np.mean(Tstat2[:,0:(j+1)],1)
        T2stats[j]=T0s[0]
        #pvalue
        pPerm0[j]= np.mean(T0s[0]<=T0s[1:numperms])
    
        #get ranks of 
        ranks=rankdata(T0s[1:numperms])
        #get empirical p-value using rank
        P0s = (numperms-ranks)/(np.float32(numperms-1))
        if j==0:
            minp0=P0s*1.0
        else:
            minp0[minp0>P0s]=P0s[minp0>P0s]*1.0
    
    pvs=max(1/np.float32(numperms),(np.sum(minp0<=np.min(pPerm0)))/np.float32(numperms-1))
    T2stats_best=T2stats[np.where(np.min(pPerm0)==pPerm0)[0][0]]
    
    
    return  pvs,T2stats_best
    
  

def sKPCA_analysis_full_data(image_names,mask,result_dir,vol_batchsize,K_components,ncore):
    
    

    images=BWAS_prepare_parallel(image_names,mask,ncore,1)
        
    
    print('Computing Laplacian matrix...')
    laplacian_mat=sKPCR_prepare_Laplacian(mask)
    laplacian_mat=sparse.csc_matrix(laplacian_mat)
        
    nsub=len(image_names)
    
    if K_components>=nsub:
        K_components=nsub-2
    
    nvol=images[0].shape[1]
        
    ind_end=int(nvol/float(vol_batchsize)-1e-4)+1
    

    for i in range(0,ind_end):
        start = time.time()
        print('Calculating the FC matrix across all subjects... '+
            str(i+1)+'/'+str(ind_end)+'...')
        
        st1=int(i*vol_batchsize)
        en1=min(nvol,int((i+1)*vol_batchsize))
        
        
        r1=np.zeros((nsub,nvol,int(en1-st1)),dtype='float32')
        for j in range(0,nsub):
            img_tmp=images[j]
            tmp=BWAS_correlation(img_tmp,img_tmp[:,st1:en1])
            tmp[tmp>0.9999]=0
            r1[j,:,:]=tmp
        
        r1=BWAS_fisher_z(r1)
        r1[np.isnan(r1)]=0
        
        r2=[np.array(a) for a in range (0,int(en1-st1))]
        for j in range(0,int(en1-st1)):
            r2[j]=r1[:,:,j];
        r1=None
        
        print('Performing sKPCA on FC matrix... '+
            str(i+1)+'/'+str(ind_end)+'...')
       
        U=np.zeros((nsub,K_components,int(en1-st1)),dtype='float32')
        
        if ncore==1:
            for j in range(0,int(en1-st1)):
                U[:,:,j]=sKPCR_linear(r2[j],laplacian_mat,K_components)                                
            
        else:                    
            num_cores = multiprocessing.cpu_count()
            Us=Parallel(n_jobs=min(num_cores-2,ncore))(delayed(
                    sKPCR_linear)(rr,laplacian_mat,K_components) for rr in r2)
            for j in range(0,int(en1-st1)):   
                U[:,:,j]=Us[j]                                
        
        print('Saving sKPCA results ...')
        np.save(result_dir+'sKPCA_map'+('%04d' % i)+'.npy',U)
        
        
        
        
        end = time.time()
        
        print('Total run time = ' + str(round(end - start,2))+ ' seconds...')
    
    
    print('Done...')
    
    return 

def sKPCR_fdr_bh(pv,thre_p):
    
    #pv = np.random.uniform(0.0, 1.0, size = (1000,))

    pv=np.sort(pv.flatten())
    
    V = len(pv);
    I = np.array(range(1,V+1),dtype='float64')

    cVID = 1
    
    tmp=np.array(np.where( pv<=(I/V*thre_p/cVID) ) )
    
    if tmp.size>0:        
        p_thre_uncorrected = pv[np.max(tmp)]
    else:
        p_thre_uncorrected=0.0

    print(p_thre_uncorrected)
    
    return p_thre_uncorrected



def sKPCR_analysis(result_dir,pheno,covariates,mask_file,vol_batchsize,numperms,ncore):
    
    
    info=nib.load(mask_file)
    mask=np.float64(info.get_data()!=0)
    affine=info.affine

    dim=np.array(np.where(mask!=0)).T
    print('Number of voxels = '+str(dim.shape[0]))
    
    nsub=pheno.shape[0]
    
    yperms = np.zeros((nsub,numperms),dtype='float32')
    yperms[:,0]=pheno.flatten()
    for j in range(1,numperms):
        yperms[:,j] = pheno[np.random.permutation(nsub)].flatten()
                
    nvol=dim.shape[0]
    
    ind_end=int(nvol/float(vol_batchsize)-1e-4)+1    


    pvs=np.ones((nvol,1))
    T2stat_best=np.zeros((nvol,1))
    for i in range(0,ind_end):
        start = time.time()
        print('Performing voxel-wise association analysis using adaptive regression... '+
            str(i+1)+'/'+str(ind_end)+'...')
        
        st1=int(i*vol_batchsize)
        en1=min(nvol,int((i+1)*vol_batchsize))
        
        pvs_tmp=np.zeros((int(en1-st1),1),dtype='float32')
        T2stat_tmp=np.zeros((int(en1-st1),1),dtype='float32')
        if ncore==1:
            U=np.load(result_dir+'sKPCA_map'+('%04d' % i)+'.npy')
            for j in range(0,int(en1-st1)):                
                pvs_tmp[j],T2stat_tmp[j]=sKPCR_adptive_regression(U[:,:,j],yperms,covariates,numperms) 
            pvs[st1:en1]=pvs_tmp 
            T2stat_best[st1:en1]=T2stat_tmp
        else:
            U=np.load(result_dir+'sKPCA_map'+('%04d' % i)+'.npy')
            num_cores = multiprocessing.cpu_count()
            PVS=Parallel(n_jobs=min(num_cores-2,ncore))(delayed(
                    sKPCR_adptive_regression)(U[:,:,j],yperms,covariates,numperms) for j in range(0,int(en1-st1)))
            for j in range(0,int(en1-st1)):   
                pvs_tmp[j]=PVS[j][0]
                T2stat_tmp[j]=PVS[j][1]
            pvs[st1:en1]=pvs_tmp 
            T2stat_best[st1:en1]=T2stat_tmp
            
        end = time.time()
        print('Number of voxels with p<=0.001 = ',str(np.sum(pvs<=0.001)))
        print('Total run time = ' + str(round(end - start,2))+ ' seconds...')
                  
      
    print('Saving result...')
    
    pval_map=np.zeros((mask.shape[0],mask.shape[1],mask.shape[2]))
    t2stat_map=np.zeros((mask.shape[0],mask.shape[1],mask.shape[2]))

    for j in range(0,nvol):
        pval_map[dim[j,0],dim[j,1],dim[j,2]]=-np.log10(pvs[j])
        t2stat_map[dim[j,0],dim[j,1],dim[j,2]]=T2stat_best[j]
    
    img = nib.Nifti1Image(pval_map, affine)
    nib.save(img, result_dir+'sKPCR_Pval_map.nii.gz')
    img = nib.Nifti1Image(t2stat_map, affine)
    nib.save(img, result_dir+'sKPCR_bestT2stat_map.nii.gz')          


    thre_p=sKPCR_fdr_bh(pvs,0.05)        
    pvs1=pvs*1.0
    pvs1[pvs1>thre_p]=1
    T2stat_best1=T2stat_best*1.0
    T2stat_best1[pvs1>thre_p]=0
    
    pval_map1=np.zeros((mask.shape[0],mask.shape[1],mask.shape[2]))
    t2stat_map1=np.zeros((mask.shape[0],mask.shape[1],mask.shape[2]))
    for j in range(0,nvol):
        pval_map1[dim[j,0],dim[j,1],dim[j,2]]=-np.log10(pvs1[j])
        t2stat_map1[dim[j,0],dim[j,1],dim[j,2]]=T2stat_best1[j]
    
    img = nib.Nifti1Image(pval_map1, affine)
    nib.save(img, result_dir+'sKPCR_Pval_map_FDR0.05.nii.gz')
    img = nib.Nifti1Image(t2stat_map1, affine)
    nib.save(img, result_dir+'sKPCR_bestT2stat_map_FDR0.05.nii.gz')


    thre_p=sKPCR_fdr_bh(pvs,0.01)        
    pvs2=pvs*1.0
    pvs2[pvs2>thre_p]=1
    T2stat_best2=T2stat_best*1.0
    T2stat_best2[pvs2>thre_p]=0
          
    pval_map2=np.zeros((mask.shape[0],mask.shape[1],mask.shape[2]))
    t2stat_map2=np.zeros((mask.shape[0],mask.shape[1],mask.shape[2]))
    for j in range(0,nvol):
        pval_map2[dim[j,0],dim[j,1],dim[j,2]]=-np.log10(pvs2[j])
        t2stat_map2[dim[j,0],dim[j,1],dim[j,2]]=T2stat_best2[j]
    
    img = nib.Nifti1Image(pval_map2, affine)
    nib.save(img, result_dir+'sKPCR_Pval_map_FDR0.01.nii.gz')
    img = nib.Nifti1Image(t2stat_map2, affine)
    nib.save(img, result_dir+'sKPCR_bestT2stat_map_FDR0.01.nii.gz')
              
    
    print('All analysis finished!')    
    
    
    
    
    return

def sKPCR_run_full_analysis(result_dir,image_dir,mask_file,toolbox_path,targets_file,cov_file,K_components,numperms,ncore):
    
    
    
    result_dir=os.path.join(result_dir,'')
    image_dir=os.path.join(image_dir,'')
    
    #get image file names
    #sort them
    img_file_names=sorted(glob.glob(image_dir+'*.nii.gz'))
    if len(img_file_names)==0:
        img_file_names=sorted(glob.glob(image_dir+'*.nii'))
        
    image_names=[]
    for i in range(0,len(img_file_names)):
        tmp=[]
        tmp.append(img_file_names[i])        
        image_names.append(tmp)
    
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
            
    nsub=targets.shape[0]
    
    #load the mask for future use
    info=nib.load(mask_file)
    mask=info.get_data()
    mask=np.float64(mask!=0)
    nvol=int(np.sum(np.sum(mask)))
    print('Number of Voxels in mask = '+str(nvol))
    print('Shape of the Mask = '+str(mask.shape))

    vol_batchsize=200
    
    #if the glm analysis is not finished yet, rerun all glm,
    #else only compute the stats
    
    #total number of voxel batchs
    ind_end=int(nvol/vol_batchsize-1e-4)+1;
    
    fname=os.path.join(result_dir,result_dir+'sKPCA_map'+('%04d' % (ind_end-1))+'.npy')
        
    if os.path.exists(fname)==0:
        sKPCA_analysis_full_data(image_names,mask,result_dir,vol_batchsize,min(K_components,nsub-5),ncore)

    sKPCR_analysis(result_dir,targets,covariates,mask_file,vol_batchsize,numperms,ncore)
    
    return










