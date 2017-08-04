function BWAS_glm(images,design,mask)
%This function performs whole-brain BWAS
%        BWAS_glm(images,design)
%Input:  images: cell (time*voxel per cell)
%        design: 2D matrix (num_sample, p variable), the first one is the
%        target variable
%Output: statistical parametric Maps saved in .mat format
%%
n_voxel=size(images{1},2);
n_sample=length(images);
n=100;
ind_end=fix(n_voxel/n)+1;
ab=dir('stat_map*.mat');
if isempty(ab)
    ind_st=1;
else
    ind_st=length(ab);
end

for i=ind_st:ind_end
    st1=(i-1)*n+1;
    en1=min(i*n,n_voxel);
    nn=n_voxel*(en1-st1+1);
    %correlation matrix...
    disp(['Calculating the correlation matrix across all subjects... '...
        num2str(i),'/',num2str(ind_end),'. '])
    r1=zeros(n_sample,n_voxel,en1-st1+1);
    for j=1:n_sample
        r1(j,:,:)=corr(images{j},images{j}(:,st1:en1));
    end
    
    r1(r1>0.9999)=0;
    rr1=0.5*log((1+r1)./(1-r1));
    rr2=reshape(rr1,[n_sample,nn]);
    
    %general linear model on FCs...
    disp(['Calculating the Statistical Parametirc Maps... '...
        num2str(i),'/',num2str(ind_end),'... '])
    ts=BWAS_Tregression(design,rr2);
    zs=norminv(tcdf(ts,n_sample-size(design,2)-1));
    stat_map=reshape(zs,[n_voxel,(en1-st1+1)]);
    
    save(['stat_map',num2str(i,'%04d'),'.mat'],'stat_map','-v7');
    
end

if exist('Estimated_fwhm.mat')==0
    disp('Estimating the smoothness of the images...')
    [d1,d2,d3]=ind2sub(size(mask),find(mask~=0));
    fwhms=[];
    for i=1:n_sample
        img=images{i};
        img=GaussianNormalization(img);
        img4d=img2dto3d(size(mask),[d1,d2,d3],img');
        fwhms(i)=mean(BWAS_est_fwhm(img4d,n_sample));
    end
    save('Estimated_fwhm.mat','fwhms');
    
end

end



