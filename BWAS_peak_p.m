function p_value1=BWAS_peak_p(mask1,mask2,t,fwhm)
%        p_value1=BWAS_peak_p(mask1,mask2,t,fwhm)


r1=(mask1/(4/3*pi))^(1/3);
R1=[1, 4*(r1^3/prod(fwhm))^(1/3), ...
    2*pi*(r1^3/prod(fwhm))^(2/3), ...
    mask1/prod(fwhm)];
r2=(mask2/(4/3*pi))^(1/3);
R2=[1, 4*(r2^3/prod(fwhm))^(1/3), ...
    2*pi*(r2^3/prod(fwhm))^(2/3), ...
    mask1/prod(fwhm)];
bonfp=1-normcdf(t)^(mask1*mask2);

EC=Gaussianfield_6Ddensity(t);

p_value=0;
for i=0:3
    for j=0:3
        p_value=p_value+R1(i+1).*R2(j+1).*EC(i+j+1,:);
    end
end

p_value1=min(p_value,1-exp(-p_value));

p_value1=min(bonfp,p_value1);


end