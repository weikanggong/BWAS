function [p1,p]=BWAS_cluster_p(mask1,mask2,CDT,clustersize,fwhm)
%        p1=BWAS_cluster_p(mask1,mask2,CDT,clustersize,fwhm)


r1=(mask1/(4/3*pi))^(1/3);
R1=[1, 4*(r1^3/prod(fwhm))^(1/3), ...
    2*pi*(r1^3/prod(fwhm))^(2/3), ...
    mask1/prod(fwhm)];
r2=(mask2/(4/3*pi))^(1/3);
R2=[1, 4*(r2^3/prod(fwhm))^(1/3), ...
    2*pi*(r2^3/prod(fwhm))^(2/3), ...
    mask1/prod(fwhm)];
V=mask1*mask2;


EC=Gaussianfield_6Ddensity(CDT);



EN=V*(normcdf(CDT,'upper'));

EL1=0;
for i=0:3
    for j=0:3
        EL1=EL1+R1(i+1).*R2(j+1).*EC(i+j+1,:);
    end
end


phi=(gamma(6/2+1)*EL1/EN)^(2/6);

p=exp(-phi.*clustersize.^(2/6));

p1=1-exp(-EL1.*p);


end