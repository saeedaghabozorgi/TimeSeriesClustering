function [details]=clustering_Hybrid_3Level(nor_traj, k,p,varargin)
plot_show=1;
options = struct('l1_dis_method','SAXminDis','l1_dtw_bound',1,'l1_rep','SAX','l1_alphabet_size',8,'l1_compression_ratio',8,'l2_dis_method','SAXminDis','l2_dtw_bound',1,'l2_rep','SAX','l2_alphabet_size',8,'l2_compression_ratio',8,'l3_dis_method','SAXminDis','l3_dtw_bound',1,'l3_rep','SAX','l3_alphabet_size',8,'l3_compression_ratio',8,'l3_alg','k-medoid');
optionNames = fieldnames(options);
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('EXAMPLE needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
    inpName = lower(pair{1}); %# make case insensitive
    if any(strmatch(inpName,optionNames))
        options.(inpName) = pair{2};
        %disp([num2str(inpName),' : ',num2str(options.(inpName))]);
        %    else
        %       error('%s is not a recognized parameter name',inpName)
    end
end
%% ----------Level 1----------------------------------------------
disp('level 1');
disp(['  ','dis_method',':',options.l1_dis_method,' ','rep',':',options.l1_rep,' ','alphabet_size',':',num2str(options.l1_alphabet_size),' ','compression_ratio',':',num2str(options.l1_compression_ratio)]);
c=[];
[c,itr]= do_kModes_time(nor_traj,k,0,'dis_method',options.l1_dis_method,'rep',options.l1_rep,'alphabet_size',options.l1_alphabet_size,'compression_ratio',options.l1_compression_ratio);
disp(['  --> DS:',num2str(length(nor_traj)),'  clusters:',num2str(k)]);
%------------------
% Evaluation
%
%   details_l1=[k,SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality];
if plot_show
    Plot_time_series_luminate(0,0,c,p,[],nor_traj,[],k,0,0.5,1);
    [SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality]= do_Evaluate(p,c,nor_traj,[],[]);
    disp(['  --> quality:',num2str(quality)]);
end;

%% ----------Level 2--CAST--------------------------------------------
disp('level 2');
disp(['  ','dis_method',':',options.l2_dis_method,' ','dtw_bound',':',num2str(options.l2_dtw_bound),' ','rep',':',options.l2_rep,' ','alphabet_size',':',num2str(options.l2_alphabet_size),' ','compression_ratio',':',num2str(options.l2_compression_ratio)]);
clusterCount=max(c);
for i=1:clusterCount
    temp_c=[];
    newData=find(c(:,1)==i);
    if length(newData)<5
        temp_c=ones(length(newData),1);
    else
        temp_c= do_CAST_time (nor_traj(newData),-1,'dis_method',options.l2_dis_method,'rep',options.l2_rep,'alphabet_size',options.l2_alphabet_size,'compression_ratio',options.l2_compression_ratio,'dtw_bound',options.l2_dtw_bound);
    end
    c(newData,2)=temp_c;
    disp(['  --> cluster:',num2str(i),'  Mems:',num2str(length(newData)),'  Clus:',num2str(max(temp_c))]);
end
%-------
% to map raw objects to new clusters
c(:,3)=c(:,1)*1000+c(:,2);
[x,y]=sort(c);
clsNum=1;
c(y(1,3),4)=clsNum;
for i=2:length(c)
    if ( x(i,3) ~= x(i-1,3))
        clsNum=clsNum+1;
        c(y(i,3),4)=clsNum;
    else
        c(y(i,3),4)=clsNum;
    end
end
l2_clusterCount=max(c(:,4));
disp(['  Number of clusters:',num2str(l2_clusterCount)]);
%--evaluation---------------------------------------------------------
if plot_show
    Plot_time_series_luminate(0,0,c(:,4),p,[],nor_traj,[],l2_clusterCount,2,0.5,2);
    % [SSEP,SSEC,RS,purity,BCubed,ConEntropy,fm]= do_Evaluate(p,cc,nor_traj,class_center,center);
    %qual_2lev=Calculate_Cluster_correct_ratio(c(:,4),p);
    %     N_reduction_2lev=1-(l2_clusterCount/rows);
    %     purity2=Calculate_Cluster_Purity(c(:,4),p,1);
end;
%% --Level 3-------------------------------------------------
% to make the prototypes
disp('level 3');
disp(['  ','dis_method',':',options.l3_dis_method,' ','dtw_bound',':',num2str(options.l3_dtw_bound),' ','rep',':',options.l3_rep,' ','alphabet_size',':',num2str(options.l3_alphabet_size),' ','compression_ratio',':',num2str(options.l3_compression_ratio)]);
disp('  Making prototype ...');
center=[];
for i=1:l2_clusterCount;
    center{i}=centre_mediod(c(:,4),i,nor_traj,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
    weight(i,1)=length(find(c(:,4)==i));
end
if plot_show 
    % Plot_time_series_luminate(0,0,c(:,4),p,center,nor_traj,[],l2_clusterCount,2,0.2,2); 
end;
disp(['  clustering:',num2str(options.l3_alg)]);
if l2_clusterCount>k
    if strmatch(options.l3_alg,'Hier_avg')
        [c3,Z]=do_Hierarchical_time(center,k,'average',-1,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
    elseif strmatch(options.l3_alg,'k-means')
        [c3,itr]= do_kMeans_time (center,k,'DTW',0,'RAW','dtw_bound',1);
    else
        [c3,~]= do_kMediod_time (center,k,0,'weight',weight, 'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
    end
    for j=1:k
        l3_mems=find(c3==j);
        for i=1:length(l3_mems)
            sub_members=find(c(:,4)==l3_mems(i));
            c(sub_members,5)=j;
        end
    end
else
    c(:,5)=c(:,4);
end

%--evaluation-------------------------------------------------------------
%if plot_show h=dendrogram(Z);end;
% Plot_time_series(0,0,c(:,5),p,[],nor_traj,t_traj,k,3,2);
%--------------------
%-- to show the centers
for j=1:length(center)
    sub_members=find(c(:,4)==j);
    xx=p(sub_members);
    pp(j,1)=mode(xx); %  to find more frequent in this cluster
end
clus_center={};
for i=1:k
    clus_center{i}=centre_mean(c(:,5),i,nor_traj);
end
[SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality]= do_Evaluate(p,c(:,5),nor_traj,[],[]);
if plot_show 
    Plot_time_series_luminate(0,0,c3,pp,[],center,[],k,2,0.5,3); 
    Plot_time_series_luminate(0,0,c(:,5),p,[],nor_traj,[],k,2,0.5,4); 
    disp(['  --> quality:',num2str(quality)]);
end;
details=[SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality];
end



function mmean=centre_mean(c,clusterNum,nor_traj)

members=find(c(:,1)==clusterNum);
n=length(nor_traj{1});
cluster_mem=zeros(size(members,1),n);
for j=1:size(members,1)
    if (members(j) > length(nor_traj))
        o=1;
    end
    cluster_mem(j,:)=nor_traj{members(j)};
end
mmean=mean(cluster_mem,1);
end


function medoid=centre_mediod(c,clusterNum,nor_traj_raw,varargin)
%clusterNum
options = struct('rep','RAW');
optionNames = fieldnames(options);
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('EXAMPLE needs propertyName/propertyValue pairs')
end
for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
    inpName = lower(pair{1}); %# make case insensitive
    if any(strmatch(inpName,optionNames))
        options.(inpName) = pair{2};
    end
end



t=find(c(:,1)==clusterNum);

if isempty (t)
    medoid=[];
elseif length(t)<=2
    nor_traj=represent_TS(nor_traj_raw(t),options.rep,varargin{:});
    % medoid=nor_traj{1}; % transfered
    medoid=nor_traj_raw{t(1)};
else
    nor_traj=represent_TS(nor_traj_raw(t),options.rep,varargin{:});
    %find distance of objects in cluster
    dis=Mtx_Distance(nor_traj,nor_traj,'same',varargin{:});
    
    %find the SSE
    dis=dis.^2;
    Error=sum(dis);
    [s,m]=min(Error);
    % medoid=nor_traj{m}; %transfered
    medoid=nor_traj_raw{t(m)};
end
end

function XXX()
% to make the prototypes
prototypes_clus=[];
center=[];
newpoints=[];
Prot_2lcluster=[];
counter=1;
% clf;
% scatter(D(:,1),D(:,2),[],'b');
minp =4;
for i=1:l2_clusterCount;
    mems=find(c(:,4)==i);
    % find rep
    if(size(mems,1)>minp)
        dis=dis_euclidean_matrix(D(mems,:),D(mems,:));
        representative=rep(dis,round(size(mems,1)/minp));
        prototypes=D(mems(representative),:);
        prototypes=shrink_tocenter(D(mems,:),prototypes,0.2);
    elseif (size(mems,1)>0)
        prototypes=mean(D(mems,:),1);
    end
    prototypes_clus{counter}=prototypes;
    newpoints=[newpoints; prototypes_clus{counter}];
    Prot_2lcluster=[Prot_2lcluster; repmat([i],size(prototypes_clus{counter},1),1)];
    counter=counter+1;
end

end

function representative=rep(dis,num)
a=zeros(size(dis,1),1);
representative=[];
U=[1:size(dis,1)]';
[m,inx]=MaxMat(dis,U);
u=U(inx);
representative(end+1,1)=u;
U(inx,:)=[];
for i=1:size(U,1)
    a(U(i),1)=a(U(i),1)+ dis(U(i),u);
end
for i=1:size(representative,1)
    a(representative(i),1)=a(representative(i),1)+ dis(u,u);
end
while size(representative,1)<num
    [~,inx]= max(a(U),[],1);
    u=U(inx);
    representative(end+1,1)=u;
    U(inx,:)=[];
    for i=1:size(U,1)
        a(U(i),1)=a(U(i),1)+ dis(U(i),u);
    end
    for i=1:size(representative,1)-1
        a(representative(i),1)=a(representative(i),1)+ dis(representative(i),u);
    end
end
end

