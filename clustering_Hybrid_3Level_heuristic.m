function [details,hurestic_param]=clustering_Hybrid_3Level_heuristic(nor_traj, k,p,dist_mtx_DTW,varargin)
disp(['-- START ---------------------------------']);
kk=round(length(nor_traj)/50);
if kk<k
    kk=k;
end

%----------------Hueristic -----------------------------
best_so_far=0;
l1_rep_arr={'SAX','RAW'};
l1_dis_method_arr={'SAXDIST','Euclid'};
l1_alg_arr={'hier_avg','k-modes','k-medoids','k-medoids-keogh'};
l3_alg_arr={'hier_avg','hier_single','k-medoids_weighted','k-medoids','k-medoids-keogh'};
for l1_rep_inx=1:length(l1_rep_arr);
    for l1_dis_method_inx=1:length(l1_dis_method_arr);
        for l1_compression_ratio=2:8
            for l1_alg_inx=1:length(l1_alg_arr);
                for l3_alg_inx=1:length(l3_alg_arr);
                    if (~(strcmp(l1_dis_method_arr{l1_dis_method_inx},'SAXDIST') && strcmp(l1_rep_arr{l1_rep_inx},'RAW'))&& ~( l1_compression_ratio>2 && strcmp(l1_rep_arr{l1_rep_inx},'RAW')) )
                        A=zeros(length(nor_traj),length(nor_traj));
                        %% ------------ Level 1-- k-mode --------
                        parameter1={'l1_dis_method',l1_dis_method_arr{l1_dis_method_inx},'l1_dtw_bound',0,'l1_rep',l1_rep_arr{l1_rep_inx},'l1_alphabet_size',8,'l1_compression_ratio',l1_compression_ratio,'l1_alg',l1_alg_arr{l1_alg_inx}};
                        [c1 ,D]=clustering_l1_preclustering(kk,p,nor_traj,0,parameter1{:});
                        %% -------------Level 2--CAST--------------------------------------------
                        parameter2={'l2_dis_method','DTW','l2_dtw_bound',1,'l2_rep','RAW','l2_alphabet_size',8,'l2_compression_ratio',2};
                        [c2 ,error_rate,N_reduction_2lev,D,A]=clustering_l2_purify(c1,p,nor_traj,D,A,0,dist_mtx_DTW,parameter2{:});
                        [center weight cen_inx]=clustering_l2_making_prototype(c2,p,nor_traj,D,0,parameter2{:});
                        %% --------------Level 3-------------------------------------------------
                        parameter3={'l3_dis_method','DTW','l3_dtw_bound',1,'l3_rep','RAW','l3_alphabet_size',8,'l3_compression_ratio',2,'l3_alg',l3_alg_arr{l3_alg_inx}};
                        [c3 ,details]=clustering_l3_merge(c2,p,k,center,cen_inx',nor_traj,weight,D,A,0,dist_mtx_DTW,parameter3{:});
                        if details(1,4)>best_so_far
                            best_so_far=details(1,4);
                            best_l1_compression_ratio=l1_compression_ratio;
                            best_l1_alg=l1_alg_arr{l1_alg_inx};
                            best_l1_rep=l1_rep_arr{l1_rep_inx};
                            best_l1_dis_method=l1_dis_method_arr{l1_dis_method_inx};
                            best_l3_alg=l3_alg_arr{l3_alg_inx};
                            
                        end
                    end
                end
            end
        end
    end
end
%---------------------------------------------------------------------
%---------------------------------------------------------------------
hurestic_param={best_so_far,best_l1_compression_ratio,best_l1_rep,best_l1_dis_method, best_l1_alg,best_l3_alg};
A=zeros(length(nor_traj),length(nor_traj));
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

