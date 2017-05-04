function  T=BWAS_peak(mask1,mask2,p_threshold,fwhm)
%         T=BWAS_peak(mask1,mask2,p_threshold,fwhm)

N=mask1*mask2;
bonfp=1-(1-p_threshold).^(1/N);
bonfT=abs(norminv(bonfp));

if p_threshold<=0.1
    T=bonfT;
else
    T=8;
end
p=BWAS_peak_p(mask1,mask2,T,fwhm)-p_threshold;
step=0.01;

while p<0
    T=T-step;
    p=BWAS_peak_p(mask1,mask2,T,fwhm)-p_threshold;
end



end