function [peak_result,cluster_result,peak_ma,cluster_ma]=BWAS_analysis_results(result_dirs,mask,peak_thre,CDT,fwhm)
%        [peak_result,cluster_result]=BWAS_analysis_results(result_dirs,mask,peak_thre,CDT,fwhm)
%Input:  result_dirs: the directory of stat_map*.mat and pca_map*.mat
%        peak_thre: between 0 and 0.1, the family-wise error rate
%                   (default 0.05, better to <0.1)
%        CDT: cluster define threshold Z
%             (default 5.5 for 3mm fmri data, better to >5)
%        fwhm: the 1*3 vector, smoothness of the image
%        n_processors: number of processors to use (default = 1)
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
    cd(result_dirs);
    peak_thre=0.05;
    CDT=5;
    load('Estimated_fwhm.mat');
    fwhm=mean(fwhms);
end

if nargin<4
    cd(result_dirs);
    CDT=5;
    load('Estimated_fwhm.mat');
    fwhm=mean(fwhms);
end
if nargin<5
    cd(result_dirs);
    load('Estimated_fwhm.mat');
    fwhm=mean(fwhms);
end

%
% if nargin<5
%     n_processors=1;
% end
% %parpool(n_processors);

if peak_thre>0.1
    disp('Peak-level inferece threshold > 0.1 is not valid...');
end

if CDT<5 && fwhm>3
    disp('Cluster-level inference is not valid with CDT(Z-statistic)<5 for FWHM>3 voxels...');
elseif CDT<4.5 && fwhm<=3
    disp('Cluster-level inference is not valid with CDT(Z-statistic)<4.5 for FWHM<3 voxels...');
end



cd(result_dirs);

a=dir('stat_map*.mat');
link_result=load(a(1).name);
[n1,n2]=size(link_result.stat_map);
%peak threshold
z_thre=BWAS_peak(sqrt(n1*(n1-1)/2),sqrt(n1*(n1-1)/2),peak_thre,[fwhm,fwhm,fwhm]+0.3);
%cluster threshold

disp('Loading GLM results...');
significant_links={};
clusters_links={};
ma1={};
parfor i=1:length(a)
    link_result=load(a(i).name);
    re=link_result.stat_map;
    re(isnan(re))=0;
    [d1,d2,v]=find(sparse(re.*(abs(re)>z_thre)));
    ma1{i}=sum(abs(re)>z_thre);
    significant_links{i}=[d1,d2+(i-1)*n2,v];
    [d1,d2,v1]=find(sparse(re.*(abs(re)>CDT)));
    clusters_links{i}=[d1,d2+(i-1)*n2,v1];
end
ma1=cell2mat(ma1);
disp('Done...');

disp('Analyzing the results...');
%% peak level inference
significant_links1=cell2mat(significant_links');
indx=significant_links1(:,1)>significant_links1(:,2);
significant_links1=significant_links1(indx,:);
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
disp(['Find ',num2str(size(peak_result,1)),' FCs survive FC-wise FWER threshold (p<',num2str(normcdf(z_thre,'upper')),')...']);


%% cluster level inference
clusters_links1=cell2mat(clusters_links');
indx=clusters_links1(:,1)>clusters_links1(:,2);
clusters_links1=clusters_links1(indx,:);
dim11=dim(clusters_links1(:,1),:);
dim22=dim(clusters_links1(:,2),:);
dis1=sqrt(sum((dim11-dim22).^2,2));
ind1=find(dis1>3);
%6d coordinates of clusters
dim111=fix(dim11(ind1,:));
dim222=fix(dim22(ind1,:));
clusters_links2=[dim111,dim222,clusters_links1(ind1,3)];
%index coordinates of clusters
clusters_links1=clusters_links1(ind1,:);

disp(['Find ',num2str(size(clusters_links2,1)),' FCs survive Cluster-defining threshold (p<',num2str(normcdf(CDT,'upper')),')...']);
%neighbourhood adjaceny matrix
adjs=sparse(pdist2(dim,dim)<=sqrt(2)+0.001);
adjs1=adjs(clusters_links1(:,1),clusters_links1(:,1));
adjs2=adjs(clusters_links1(:,2),clusters_links1(:,2));
adjs3=(adjs1==1 & adjs2==1);

%connected components
[comps,compsize]=get_components(adjs3);

disp(['Find ',num2str(length(compsize)),' FC clusters...']);

cluster_result_tmp=[];
indss={};

for jj=1:length(compsize)
    indx=(comps==jj);
    indss{jj}=indx;
    voxel_1d_indx=clusters_links1(indx,:);
    cluster_result_tmp(jj,1)=compsize(jj);
    [cluster_result_tmp(jj,3),cluster_result_tmp(jj,2)]=BWAS_cluster_p(sqrt(n1*(n1-1)/2),sqrt(n1*(n1-1)/2),CDT,compsize(jj),[fwhm,fwhm,fwhm]+0.3);
    cluster_result_tmp(jj,4)=length(unique(voxel_1d_indx(:,1)));
    cluster_result_tmp(jj,5)=length(unique(voxel_1d_indx(:,2)));
    
end

cluster_result_tmp(:,6)=mafdr(cluster_result_tmp(:,2),'BHFDR',1);
[cluster_result_tmp,ind3]=sortrows(cluster_result_tmp,2);

disp([num2str(sum(cluster_result_tmp(:,3)<0.05)),' FC clusters survive FWER-corrected p<0.05 under this CDT...']);

indss=indss(ind3);
cluster_result=struct;
for jj=1:length(compsize)
    cluster_result(jj).clusters=clusters_links2(indss{jj},:);
    cluster_result(jj).cluster_size=cluster_result_tmp(jj,1);
    cluster_result(jj).uncorrected_p=cluster_result_tmp(jj,2);
    cluster_result(jj).FWER_p=cluster_result_tmp(jj,3);
    cluster_result(jj).region1_nvoxels=cluster_result_tmp(jj,4);
    cluster_result(jj).region2_nvoxels=cluster_result_tmp(jj,5);
    cluster_result(jj).cluster_FDR=cluster_result_tmp(jj,6);
end

% img=img2dto3dmask([61,73,61],cluster_result(2).clusters(:,1:3));
% cc=bwconncomp(img,18);
% cc

%mapping peaks and clusters to brain region
labels=unique(mask(:));
labels(labels==0)=[];
if length(labels)>1
    %peak mapping
    kk1=size(peak_result,1);
    if kk1>30
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

disp('The analysis is finished!');

save(['Link_BWAS_results_CDT',num2str(CDT),'.mat'],'peak_result','cluster_result','mask','peak_thre','CDT','fwhm','peak_ma','cluster_ma')


return


function [comps,comp_sizes] = get_components(adj)
%GET_COMPONENTS     connected components
%
%   [comps,comp_sizes] = get_components(adj);
%
%   Returns the components of an undirected graph specified by the binary and
%   undirected adjacency matrix adj. Components and their constitutent nodes are
%   assigned the same index and stored in the vector, comps. The vector, comp_sizes,
%   contains the number of nodes beloning to each component.
%
%   Inputs:         adj,    binary and undirected adjacency matrix
%
%   Outputs:      comps,    vector of component assignments for each node
%            comp_sizes,    vector of component sizes
%
%   Note: disconnected nodes will appear as components with a component
%   size of 1
%
%   J Goni, University of Navarra and Indiana University, 2009/2011


% if size(adj,1)~=size(adj,2)
%     error('this adjacency matrix is not square');
% end
%
% if ~any(adj-triu(adj))
%     adj = adj | adj';
% end
%
% %if main diagonal of adj do not contain all ones, i.e. autoloops
% if sum(diag(adj))~=size(adj,1)
%
%     %the main diagonal is set to ones
%     adj = adj|speye(size(adj));
% end

%Dulmage-Mendelsohn decomposition
[~,p,~,r] = dmperm(adj);

%p indicates a permutation (along rows and columns)
%r is a vector indicating the component boundaries

% List including the number of nodes of each component. ith entry is r(i+1)-r(i)
comp_sizes = diff(r);

% Number of components found.
num_comps = numel(comp_sizes);

% initialization
comps = zeros(1,size(adj,1));

% first position of each component is set to one
comps(r(1:num_comps)) = ones(1,num_comps);

% cumulative sum produces a label for each component (in a consecutive way)
comps = cumsum(comps);

%re-order component labels according to adj.
comps(p) = comps;

return










