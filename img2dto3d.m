function img=img2dto3d(size3dimage,dim_vector,value)
 m=size(value,2);
 img=zeros([size3dimage,m]);
 n=size(dim_vector,1);
 for i=1:n
     img(dim_vector(i,1),dim_vector(i,2),dim_vector(i,3),:)=value(i,:);
 end
end

