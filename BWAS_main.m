function BWAS_main(result_dir,image_dir,design_dir,mask_dir,CDT,FWER_p)
%        BWAS_main(result_dir,image_dir,design_dir,mask_dir,CDT,cluster_p)

if nargin<5
    CDT=5;
    FWER_p=0.05;
end


%% prepare the data

if exist('data_to_use.mat')==0
    cd(image_dir);
    image_dir1=dir('*.nii.gz');
    names={};
    for i=1:length(image_dir1)      
        names{i}=image_dir1(i).name;
    end
   
    aa=load_nii(mask_dir);
    mask=aa.img;
    images=BWAS_prepare(names,mask);

    disp('Checking image quality...');
    n_voxel=size(images{1},2);
    n_sample=length(images);
    bad_voxels={};
    n_bad=0;
    for i=1:n_sample
        images{i}(isnan(images{i}))=0;
        ind=all(images{i})<0.0001;
        bad_voxels{i}=find(ind==1);
        n_bad=n_bad+sum(ind);
        images{i}(:,ind)=0;
    end
    if n_bad>0
        disp(['Your images contrains ',num2str(n_bad),...
            ' "bad" voxels, please check "bad_voxels.mat"!'])
        save('bad_voxels.mat','bad_voxels');
    else
        disp('No "bad" voxels are found.')
    end
    
    design=readtable(design_dir);
    design=table2array(design);
    
    disp('Checking phenotype quality...');
    if class(design)~='double'
        error('The data type of design matrix should be double!')
    end
    if size(design,1)~=n_sample
        error('The number of rows in the design matrix is not the same as the number of images!')
    end
    if sum(isnan(design(:)))>0
        error('Your phenotype matrix contains NaN, please check it!');
    end
    disp('finished! ')
    
    bad_voxels1=unique(cell2mat(bad_voxels));
    [d1,d2,d3]=ind2sub(size(mask),find(mask~=0));
    dimdim=[d1,d2,d3];
    mask_new=img2dto3dmask(size(mask),dimdim(setdiff(1:n_voxel,bad_voxels1),:));
    mask=mask_new;
    for i=1:n_sample
        images{i}=images{i}(:,setdiff(1:n_voxel,bad_voxels1));
    end
    disp('Saving the data for future analysis...')
    save('data_to_use.mat','images','design','mask','-v7.3')
    
else
    cd(image_dir);
    disp('Loading the data...');
    load('data_to_use.mat');
end

%% perform the analysis
cd(result_dir);
mkdir results
cd results
BWAS_glm(images,design,mask);
cd ..
%% analysis the results
[peak_result,cluster_result,peak_ma,cluster_ma]=BWAS_analysis_link('results',mask,FWER_p,CDT);

aa=load_nii(mask_dir);
aa.img=peak_ma;
save_nii(aa,'peak_MA.nii.gz')
aa.img=cluster_ma;
save_nii(aa,'cluster_MA.nii.gz')
movefile('Link_BWAS_results.mat','..');
movefile('peak_MA.nii.gz','..');
movefile('cluster_MA.nii.gz','..');

end
