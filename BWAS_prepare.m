function images=BWAS_prepare(imgs_abs_dir,mask)
%        BWAS_prepare(imgs_abs_dir,mask,name)
%Input   imgs_abs_dir: cell array, one cell one absolute path of preprocessed fmri image
%        mask: 3D matrix, non-zero elements are the voxel you want to use
%        in the analysis.
%Output  images: the data can be used in subsequent analysis.

images={};
parfor i=1:length(imgs_abs_dir)
    info=load_nii(imgs_abs_dir{i});
    img=info.img;
    img2=[];
    for j=1:size(img,4)
        img1=img(:,:,:,j);
        img2(j,:)=img1(mask~=0);
    end
    images{i}=img2;
    i
end



end
