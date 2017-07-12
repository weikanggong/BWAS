function BWAS_glm(images,design,mask)
%This function performs whole-brain BWAS
%        BWAS_glm(images,design,mask)
%Input:  images: cell (time*voxel per cell)
%        design: 2D matrix (num_sample, p variable), the first one is the
%        target variable
%        mask: the voxel you want to analyse is 1 , otherwise 0.
%Output: Statistical parametric Maps saved in .mat format
%%
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
    disp(['Your images contain ',num2str(n_bad),...
        ' "bad" voxels, please check "bad_voxels.mat"!'])
    save('bad_voxels.mat','bad_voxels');
else
    disp('No "bad" voxels are found.')
end

disp('Checking phenotype quality...');
if class(design)~='double'
    error('The data type of design matrix should be double!')
end
if sum(isnan(design(:)))>0
    error('Your phenotype matrix contains NaN, please check it!');
end
disp('finished! ')

n=100;
ind_end=fix(n_voxel/n)+1;


tic;
for i=1:ind_end
    st1=(i-1)*n+1;
    en1=min(i*n,n_voxel);
    nn=n_voxel*(en1-st1+1);
    %correlation matrix...
    disp(['Calculating the correlation matrix across all samples... '...
        num2str(i),'/',num2str(ind_end),'. '])
    r1=zeros(n_sample,n_voxel,en1-st1+1);
    for j=1:n_sample
        r1(j,:,:)=corr(images{j},images{j}(:,st1:en1));
    end

    r1(r1>0.9999)=0;
    rr1=0.5*log((1+r1)./(1-r1));
    rr2=reshape(rr1,[n_sample,nn]);
    
    %general linear model on FCs...
    disp(['Calculating the Statistical Parametric Maps... '...
        num2str(i),'/',num2str(ind_end),'. '])
    ts=BWAS_Tregression(design,rr2);
    stat_map=reshape(ts,[n_voxel,(en1-st1+1)]);
    
    save(['stat_map',num2str(i,'%04d'),'.mat'],'stat_map','-v7');
       
end
tt=toc;

disp('Estimating smoothness of the images...')
[d1,d2,d3]=ind2sub(size(mask),find(mask~=0));
fwhms=[];
for i=1:n_sample
    img=images{i};
    img=GaussianNormalization(img);
    img4d=img2dto3d(size(mask),[d1,d2,d3],img');
    fwhms(i)=mean(BWAS_est_fwhm(img4d,n_sample));
end
save('Estimated_fwhm.mat','fwhms');

disp(['Done... Total running time = ',num2str(tt),' seconds']);
end



