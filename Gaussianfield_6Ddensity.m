function EC=Gaussianfield_6Ddensity(t)

EC=[];

a = 4.*log(2);
a1= 2.*pi;
b = exp(-t.^2/2);

EC(1,:) = normcdf(t,'upper');
EC(2,:) = a.^(1/2)/a1.*b;
EC(3,:) = a./(a1.^(3/2)).*b.*t;
EC(4,:) = a.^(3/2)/(a1.^2).*b.*(t.^2 - 1);
EC(5,:) = a.^2/(a1.^(5/2)).*b.*(t.^3-3.*t);
EC(6,:) = a.^(5/2)/(a1.^3).*b.*(t.^4-6.*t.^2+3);
EC(7,:) = a.^3/(a1.^(7/2)).*b.*(t.^5-10.*t.^3+15.*t);

end