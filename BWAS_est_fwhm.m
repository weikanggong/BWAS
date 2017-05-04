function fwhm = BWAS_est_fwhm(X,df)
%
%        fwhm = BWAS_est_fwhm(X,df)
% input:
%    X:        A 3(or 4)-Dimensional matrix of statistic image (the 4-th dimension contains multiple images).
%    df:       the degree of freedom of general linear model
% output:
%    fwhm:  the FWHM parameter in x,y,z dimension
% details:
%    The fwhm value is derived from the roughness matrix Lambda. Lambda is the
%    correlation matrix of gradient of a Gaussian field in x, y, and z
%    directions, and in a typical SPM analysis Lambda is derived from
%    residual images. In this function, Lambda is derived from a statistic
%    image directly using the theoretical expression of the grandient of a
%    T-field and an F-field in [1]. This can be done by scaling the covariance
%    matrix of numerical grandients of a statistic image appropriately. Based
%    on Lambda, fwhm is calculated and returned as an output of this function.
%
%
% REFERENCE:
% [1]. Worsley KJ. 
%         Local maxima and the expected Euler characteristic of excursion sets
%         of chi-square, F and t fields.
%         Advances in Applied Probability, 26: 13-42 (1994)
%
% Modified from pm_est_fwhm
[n1,n2,n3,n4]=size(X);

muX   = df^(1/2)*(df-1)/(df-2);
varX  = 2*df*(df-1)/((df-2)^2*(df-4));
muY      = 2^(-1/2)*exp(gammaln((df-1)/2) - gammaln(df/2));
varY     = 1/(df-2) - muY^2;

%-scaling factor for var(derivative) matrix
Dscale   = 1/(varX*varY + muX^2*varY + muY^2*varX + muX^2*muY^2);

%-NaN masking the image
X      = X.*X./X;

%-Deriv in x direction
dx = diff(X,1,1); %-Deriv in x direction
varX=var(reshape(dx,[(n1-1)*n2*n3,n4]),'omitnan');

%-Deriv in y direction
dy = diff(X,1,2); %-Deriv in y direction
varY=var(reshape(dy,[n1*(n2-1)*n3,n4]),'omitnan');

%-Deriv in z direction
dz = diff(X,1,3); %-Deriv in z direction
varZ=var(reshape(dz,[n1*n2*(n3-1),n4]),'omitnan');

%-calculating global FWHM
fwhm     = (4*log(2))^(1/2)*([varX.*varY.*varZ].*(Dscale^3)).^(-1/6);

			    

