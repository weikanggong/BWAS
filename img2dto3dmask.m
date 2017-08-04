function img=img2dto3dmask(size3dimage,dim_vector)
 img=zeros(size3dimage);
 for i=1:size(dim_vector,1)
     img(dim_vector(i,1),dim_vector(i,2),dim_vector(i,3))=1;
 end
end