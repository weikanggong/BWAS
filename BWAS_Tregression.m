%      T=BWAS_Tregression(X,Y)
%Input:
%      X: n-by-p matrix, n sample and p covariates (The first column is the target variable)
%      Y: n-by-m matrix, n sample and m voxels
%Output: T statistics and Residuals
function [T,Res]=BWAS_Tregression(X,Y)

    contrast=[1,zeros(1,size(X,2))]';
    
    X=[X,ones(size(X,1),1)];

    df=size(X,1)-size(X,2);
    
    ss=inv(X'*X);
    
    beta=ss*X'*Y;
    
    Res=Y-X*beta;
      
    sigma=sqrt(sum(Res.^2)./df);
    
    T=beta'*contrast./(sigma'*sqrt(contrast'*ss*contrast));
     
    %stdRes=Res./sigma;
    
    %p=2*(1-tcdf(abs(T),df));
    


end


