function [details]=clustering_Hybrid_3Level(nor_traj, k,p,varargin)
plot_show=0;
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

  %% ------------ Level 1-- k-mode --------
    parameter1={'l1_dis_method','Euclid','l1_dtw_bound',1,'l1_rep','RAW','l1_alphabet_size',8,'l1_compression_ratio',6};
    c=clustering_l1_preclustering(k,p,nor_traj,parameter1{:});
   
    %% -------------Level 2--CAST--------------------------------------------
    parameter2={'l2_dis_method','DTW','l2_dtw_bound',0.8,'l2_rep','RAW','l2_alphabet_size',8,'l2_compression_ratio',2};
    c=clustering_l2_purify(c,p,nor_traj,parameter2{:});
    [center weight]=clustering_l2_making_prototype(c,nor_traj,parameter2{:});
    %% --------------Level 3-------------------------------------------------
    parameter3={'l3_dis_method','DTW','l3_dtw_bound',0.9,'l3_rep','RAW','l3_alphabet_size',8,'l3_compression_ratio',6,'l3_alg','k-medoid'};
    [c ,details]=clustering_l3_merge(c,p,k,center,nor_traj,weight,parameter3{:});

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

