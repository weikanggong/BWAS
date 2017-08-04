function [peak_result,cluster_result,peak_ma,cluster_ma]=BWAS_analysis_link(result_dirs,mask,peak_thre,CDT)

%        results=BWAS_analysis_link(result_dirs,peak_thre,CDT)
%Input:  result_dirs: the directory of stat_map*.mat and pca_map*.mat
%        peak_thre: between 0 and 0.1, the family-wise error rate
%                   (default 0.05, better to <0.1)
%        CDT: cluster define threshold Z
%             (default 5.5 for 3mm fmri data, better to >5)
%        fwhm: the 1*3 vector, smoothness of the image
%Output: save as 'Link_BWAS_results.mat' in 'result_dirs':
%        peak_result (table): The results for peak-level inference
%                             Each row is a significant functional
%                             connectivity
%                             The first six columns if the matrix
%                             coordinates of functional connectivities,
%                             The seventh column is the Z-statistics of the
%                             functional connectivities (two-sided).
%        cluster_result(structure): The results for cluster-level
%                                   inference. Each structure is a
%                                   functional-connectivity cluster.
%             cluster_results.clusters: the functional connectivities within the clusters.
%                                       The first six columns if the matrix
%                                       coordinates of functional connectivities,
%                                       The seventh column is the Z-statistics of the
%                                       functional connectivities (two-sided).
%             cluster_results.cluster_size: the number of functional connectivities in this cluster
%             cluster_result.uncorrected_p: the uncorrected p-value of this
%                                           cluster.
%             cluster_result.FWER_p; the FWER corrected p-value of this
%                                        cluster.
%             cluster_result.region1_nvoxels: The number of voxels in the
%                                             first endpoint-cluster of this functional connectivity
%                                             cluster.
%             cluster_result.region2_nvoxels:The number of voxels in the
%                                             second endpoint-cluster of this functional connectivity
%                                             cluster.
%        fwhm: The estimated mean FWHM of your data
%        peak_ma: A 3D matrix (the size is the same as your mask), each voxel
%                 is the number of significant functional connectivities for
%                 peak-level inference.
%        cluster_ma: A 3D matrix (the size is the same as your mask), each voxel
%                 is the number of significant functional connectivities for
%                 cluster-level inference.


%% begin analysis
if nargin<3
    peak_thre=0.05;
    CDT=5;
end

if nargin<4
    CDT=5;
end


if peak_thre>0.1
    disp('Not recommanded for peak threshold > 0.01...');
end

if CDT<5
    disp('Not recommanded for CDT<5...');
end



cd(result_dirs);
load('Estimated_fwhm.mat');
fwhm=mean(fwhms);

a=dir('stat_map*.mat');
link_result=load(a(1).name);
[n1,n2]=size(link_result.stat_map);
%peak threshold
z_thre=BWAS_peak(sqrt(n1*(n1-1)/2),sqrt(n1*(n1-1)/2),peak_thre,[fwhm,fwhm,fwhm]);
%cluster threshold

disp('Loading the SPMs...');
significant_links={};
clusters_links={};
ma1={};ma2={};
parfor i=1:length(a)
    link_result=load(a(i).name);
    re=link_result.stat_map;
    re(isnan(re))=0;
    [d1,d2,v]=find(sparse(re.*(abs(re)>z_thre)));
    ma1{i}=sum(abs(re)>z_thre);
    significant_links{i}=[d1,d2+(i-1)*n2,v];
    [d1,d2,v1]=find(sparse(re.*(abs(re)>CDT)));
    clusters_links{i}=[d1,d2+(i-1)*n2,v1];
    ma2{i}=sum(abs(re)>CDT);
end
ma1=cell2mat(ma1);
ma2=cell2mat(ma2);


significant_links1=cell2mat(significant_links');
indx=significant_links1(:,1)>significant_links1(:,2);
significant_links1=significant_links1(indx,:);

clusters_links1=cell2mat(clusters_links');

disp('Analyzing the results...');
%% peak level inference
[d1,d2,d3]=ind2sub(size(mask),find(mask~=0));
dim=[d1,d2,d3];
dim1=dim(significant_links1(:,1),:);
dim2=dim(significant_links1(:,2),:);
%discard nearby FCs
dis=sqrt(sum((dim1-dim2).^2,2));
ind=find(dis>5);

significant_links2=[fix(dim1(ind,:)),fix(dim2(ind,:)),significant_links1(ind,3)];

%peak level inference result
peak_result=array2table(significant_links2);
peak_result.Properties.VariableNames={'x1','x2','x3','x4','x5','x6','Z_statistics'};
disp(['Find ',num2str(size(peak_result,1)),' significant peak FCs...']);

%% cluster level inference
dim11=dim(clusters_links1(:,1),:);
dim22=dim(clusters_links1(:,2),:);
dis1=sqrt(sum((dim11-dim22).^2,2));
ind1=find(dis1>5);
dim111=fix(dim11(ind1,:));
dim222=fix(dim22(ind1,:));


clusters_links2=[dim111,dim222,clusters_links1(ind1,3)];

%find their mask
dim1111=sub2ind(size(mask),dim111(:,1),dim111(:,2),dim111(:,3));
dim2222=sub2ind(size(mask),dim222(:,1),dim222(:,2),dim222(:,3));
mask11=img2dto3dmask(size(mask),clusters_links2(:,1:3));
mask22=img2dto3dmask(size(mask),clusters_links2(:,4:6));
mask33=img2dto3d(size(mask),dim,ma2');
mask11=(mask11 & mask33>3);
mask22=(mask22 & mask33>3);
%find connected components
cc1=bwconncomp(mask11,18);
cc2=bwconncomp(mask22,18);
l1=cc1.PixelIdxList;
l2=cc2.PixelIdxList;

%find number of links between each connected components
ind1=[];
ind2=[];
parfor j=1:size(dim1111,1)
    index1 = cellfun(@(x) max(x==dim1111(j)), l1, 'UniformOutput', 1);
    index2 = cellfun(@(x) max(x==dim2222(j)), l2, 'UniformOutput', 1);
    if sum(index1)~=0 && sum(index2)~=0
        ind1(j)=find(index1==1);
        ind2(j)=find(index2==1);
    end
end

n_links=zeros(length(l1),length(l2));
for k1=1:length(l1)
    for k2=k1:length(l2)
        n_links(k1,k2)=sum(ind1==k1 & ind2==k2);
    end
end

%the number of links (clsize) and their clusters
[dd1,dd2,clsize]=find(sparse(n_links));
[clsize,ind11]=sort(clsize,'descend');
dd1=dd1(ind11);dd2=dd2(ind11);
disp(['Find ',num2str(size(dd1,1)),' FC clusters...']);

cluster_result=struct;
rawp=[];
for jj=1:length(dd1)
    ind3=(ind1==dd1(jj) & ind2==dd2(jj));
    clsize1=clsize(jj);
    %cluster uncorrected and FWER-corrected p-value
    [p1,rawp(jj)]=BWAS_cluster_p(sqrt(n1*(n1-1)/2),sqrt(n1*(n1-1)/2),CDT,clsize1,[fwhm,fwhm,fwhm]);
    if clsize1>0
        %FC in the cluster
        cluster_result(jj).clusters=clusters_links2(ind3,:);
        %cluster size
        cluster_result(jj).cluster_size=clsize1;
        
        cluster_result(jj).uncorrected_p=rawp(jj);
        cluster_result(jj).FWER_p=p1;
        cluster_result(jj).region1_nvoxels=length(l1{dd1(jj)});
        cluster_result(jj).region2_nvoxels=length(l2{dd2(jj)});
    end
end
%cluster's false discovery rate
if ~isempty(rawp)
    fdr=mafdr(rawp,'BHFDR',1);
    for jj=1:length(dd1)
        if clsize(jj)>0
            cluster_result(jj).FDR=fdr(jj);
        end
    end
end
%mapping peaks and clusters to brain region
labels=unique(mask(:));
labels(labels==0)=[];
if length(labels)>1
    %peak mapping
    kk1=size(peak_result,1);
    if kk1>0
        region1=[];
        region2=[];
        for i=1:kk1
            region1(i)=mask(peak_result.x1(i),peak_result.x2(i),peak_result.x3(i));
            region2(i)=mask(peak_result.x4(i),peak_result.x5(i),peak_result.x6(i));
        end
        peak_result.region1=region1';
        peak_result.region2=region2';
    end
    %cluster mapping
    kk2=size(cluster_result,2);
    for i=1:kk2
        region1=[];
        region2=[];
        dat=cluster_result(i).clusters;
        for j=1:size(dat,1)
            region1(j)=mask(dat(j,1),dat(j,2),dat(j,3));
            region2(j)=mask(dat(j,4),dat(j,5),dat(j,6));
        end
        una=unique(region1);
        prop=histc(region1,una)/length(region1);
        [prop1,ind]=sort(prop,'descend');
        una1=una(ind);
        k=find(cumsum(prop1)>0.9,1);
        xx={};
        for ii=1:k
            xx{ii,1}=[num2str(prop1(ii)*100),'% FCs connect region ',num2str(una1(ii))];
        end
        
        una=unique(region2);
        prop=histc(region2,una)/length(region2);
        [prop1,ind]=sort(prop,'descend');
        una1=una(ind);
        k=find(cumsum(prop1)>0.9,1);
        for ii=1:k
            xx{ii,2}=[num2str(prop1(ii)*100),'% FCs connect region ',num2str(una1(ii))];
        end
        xx1=cell2table(xx);
        xx1.Properties.VariableNames={'FC_start_region','FC_end_region'};
        cluster_result(i).region_label=xx1;
    end
    
end

peak_ma=img2dto3d(size(mask),[d1,d2,d3],ma1');
cluster_ma=zeros(size(mask));
if ~isempty(rawp)
    for j=1:length(cluster_result)
        if cluster_result(j).FWER_p<0.05
            clre=cluster_result(j).clusters;
            d1=clre(:,1);
            d2=clre(:,2);
            d3=clre(:,3);
            d4=clre(:,4);
            d5=clre(:,5);
            d6=clre(:,6);
            for i=1:length(d1)
                cluster_ma(d4(i),d5(i),d6(i))=cluster_ma(d4(i),d5(i),d6(i))+1;
                cluster_ma(d1(i),d2(i),d3(i))=cluster_ma(d1(i),d2(i),d3(i))+1;
            end
        end
    end
end
disp('The analysis is finished!');

save('Link_BWAS_results.mat','peak_result','cluster_result','mask','peak_thre','CDT','fwhm','peak_ma','cluster_ma')


end