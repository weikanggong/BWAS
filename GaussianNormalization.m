%data: rows are features, columns are samples
%normalize each column to mean 0 variance 1.

function new_data=GaussianNormalization(data)

m1=mean(data);

v1=var(data);

n=size(data,1);

new_data=(data-repmat(m1,n,1))./repmat(sqrt(v1),n,1);



end