function [c,itr]= do_kMeans_time (nor_traj,k,dis_method,isRand,rep,varargin)

options = struct('alphabet_size',0,'compression_ratio',0);
optionNames = fieldnames(options);
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('EXAMPLE needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
   inpName = lower(pair{1}); %# make case insensitive
   if any(strmatch(inpName,optionNames))
      options.(inpName) = pair{2};
%    else
%       error('%s is not a recognized parameter name',inpName)
    end
end


c=[];
Rows=length(nor_traj);
temp=zeros(Rows,1);
itr=0;

if nargin<3,        isRand=0;   end
if nargin<2,        k=1;        end
if nargin<1 end


if k>Rows
    return;
end
% initial value of centroid
if isRand,
    p = randperm(Rows);      % random initialization
    for i=1:k
        center{i}=nor_traj{p(i)};
    end
else
    for i=1:k
        center{i}=nor_traj{i};      % sequential initialization
    end
end

while 1,
    itr=itr+1;
    itr
    dis=Mtx_Distance(nor_traj,center,'cell_not_same',dis_method,varargin{:});
    [z,c]=min(dis,[],2);  % find group matrix g
    if (c==temp | itr==10),
        break;          % stop the iteration
    else
        temp=c;         % copy group matrix to temporary variable
    end
    
    switch rep
        case 'SAX'
            %-------------------------------------------------------------
            for i=1:k
                 center{i}= centre_modes(c,i, nor_traj);
                % center{i}= medoid_SAX(c,i, nor_traj,dis_method,alphabet_size,compression_ratio);
               % center{i}=centroid_SAX(c,i,nor_traj,length(nor_traj{1}), alphabet_size,compression_ratio)
                if  isempty(center{i})
                    center{i}=nor_traj{randperm(Rows)};
                end
            end
            %-------------------------------------------------------------
        case 'RAW'
            for i=1:k
                center{i}=centre_mean(c,i,nor_traj);
            end
            %-------------------------------------------------------------
        case 'PAA'
            for i=1:k
                center{i}=centre_mean(c,i,nor_traj);
            end
            %-------------------------------------------------------------
        otherwise
            error(sprintf('DOCLUSTERING - unsupported algorithm "%s"',dis_method))
    end %of switch/case
end
end

function mmean=centre_mean(c,clusterNum,nor_traj)
t=find(c(:,1)==clusterNum);
if isempty(t)
    p = randperm(length(nor_traj));      % random initialization
    mmean= nor_traj{p(1)};
else
    n=length(nor_traj{1});
    cluster_mem=zeros(1,n);
    for j=1:size(t,1)
        cluster_mem(j,:)=nor_traj{t(j)};
    end
    mmean=mean(cluster_mem,1);
end
end

function medoid=medoid_SAX(c,clusterNum,SAX_nor_traj,dis_method,alphabet_size,compression_ratio)
t=find(c(:,1)==clusterNum);
if isempty(t)
    mmean=[];
else
    dis=Mtx_Distance(SAX_nor_traj(t),SAX_nor_traj(t),'same',dis_method,alphabet_size,compression_ratio);
    dis=dis^2;
    sum_dis=sum(dis);
    [s,m]=min(sum_dis);
    mmean=SAX_nor_traj{t(m)};
end
medoid=mmean;
end

function centroid=centroid_SAX(c,clusterNum,nor_traj,data_len, alphabet_size,compression_ratio)
t=find(c(:,1)==clusterNum);
if isempty(t)
    p = randperm(length(nor_traj));      % random initialization
    mmean= nor_traj{p(1)};
else
    n=length(nor_traj{1});
    cluster_mem=zeros(1,n);
    for j=1:size(t,1)
        cluster_mem(j,:)=nor_traj{t(j)};
    end
    mmean=mean(cluster_mem,1);
end
nseg=data_len/compression_ratio;
centroid= rep_SAX(mmean, data_len, nseg, alphabet_size);
end

function modes=centre_modes(c,clusterNum,nor_traj)
t=find(c(:,1)==clusterNum);
if isempty (t)
    modes=[];
elseif length(t)==1
    modes=nor_traj{t(1)};
else
   
    for i=1:length(t)
        X(i,:)= nor_traj{t(i)};
    end
    modes=mode(X);
end
end
