function fwhm = BWAS_est_fwhm(X)
%
%        fwhm = BWAS_est_fwhm(X,df)
% Input:
%    X:   A 3(or 4)-Dimensional matrix of statistic image (the 4-th dimension are times).
% Output: Estimated fwhm of each image.

[n1,n2,n3,n4]=size(X);


X = X.*X./X;

dx = diff(X,1,1); 
varX=var(reshape(dx,[(n1-1)*n2*n3,n4]),'omitnan');

dy = diff(X,1,2); 
varY=var(reshape(dy,[n1*(n2-1)*n3,n4]),'omitnan');

dz = diff(X,1,3); 
varZ=var(reshape(dz,[n1*n2*(n3-1),n4]),'omitnan');

varXYZ=[varX.*varY.*varZ].^(1/3);

varImg=[];
for i=1:n4
    tmp=X(:,:,:,i);
    varImg(i)=var(tmp(~isnan(tmp)));
end

fwhm=real(sqrt(-2*log(2)./log(1-varXYZ/2./varImg)));
			   
