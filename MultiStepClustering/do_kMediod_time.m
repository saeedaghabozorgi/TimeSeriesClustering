function [c,itr]= do_kMediod_time (nor_traj_raw,k,isRand,dist_mtx_DTW,varargin)

options = struct('alphabet_size',0,'compression_ratio',0,'weight',[],'rep','RAW');
optionNames = fieldnames(options);
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('EXAMPLE needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
    inpName = lower(pair{1}); %# make case insensitive
    if any(strmatch(inpName,optionNames))
        options.(inpName) = pair{2};
        % disp(['  ', num2str(inpName),' : ',num2str(options.(inpName))]);
        %    else
        %       error('%s is not a recognized parameter name',inpName)
    end
end


% representation
nor_traj=represent_TS(nor_traj_raw,options.rep,varargin{:});

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
    inx = randperm(Rows);      % random initialization
    center=nor_traj{inx(1:k)};
else
    inx=[1:k];
    center=nor_traj(inx);      % sequential initialization
end

if ~isempty(options.weight)
    [~,inx]=sort(options.weight,'descend');
    inx=inx';
    center=nor_traj(inx(1:k));
end
cen_inx=inx(1:k);

while 1,
    itr=itr+1;
    if isempty(dist_mtx_DTW)
        dis=Mtx_Distance(nor_traj,center,'cell_not_same','Norm',varargin{:});
    else
        dis= dist_mtx_DTW(:,cen_inx);
    end
    [z,c]=min(dis,[],2);  % find group matrix g
    if (c==temp | itr==20),
        break;          % stop the iteration
    else
        temp=c;         % copy group matrix to temporary variable
    end
    
    switch options.rep
        %-----------------------------------------------------------------
        case 'SAX'
            for i=1:k
                if isempty(options.weight)
                      newData=find(c(:,1)==i);
                    dis1=dist_mtx_DTW(newData,newData);
                    [medoid inx]=centre_mediod(c,i,nor_traj,dis1,varargin{:});
                    cen_inx(1,i)=inx;
                    center{i}=medoid;
                else
                        newData=find(c(:,1)==i);
                    dis1=dist_mtx_DTW(newData,newData);
                    [medoid inx]=centre_mediod_weighted(c,i,nor_traj,dis1,varargin{:});
                    cen_inx(1,i)=inx;
                    center{i}=medoid;
                end
                % center{i}=centroid_SAX(c,i,nor_traj,length(nor_traj{1}), alphabet_size,compression_ratio)
                if  isempty(center{i})
                    center{i}=nor_traj{randperm(Rows)};
                end
            end
            %-------------------------------------------------------------
        case 'RAW'
            for i=1:k
                if isempty(options.weight)
                    newData=find(c(:,1)==i);
                    dis1=dist_mtx_DTW(newData,newData);
                    [medoid inx]=centre_mediod(c,i,nor_traj,dis1,varargin{:});
                    cen_inx(1,i)=inx;
                    center{i}=medoid;
                else
                    newData=find(c(:,1)==i);
                    dis1=dist_mtx_DTW(newData,newData);
                    [medoid inx]=centre_mediod_weighted(c,i,nor_traj,dis1,varargin{:});
                    cen_inx(1,i)=inx;
                    center{i}=medoid;
                    
                end
                
                
                
            end
            %-------------------------------------------------------------
        case 'PAA'
            for i=1:k
                [medoid inx]=centre_mediod(c,i,nor_traj,dis,varargin{:});
                cen_inx(i,1)=inx;
                center{i}=medoid;
            end
            %-------------------------------------------------------------
        otherwise
            error(sprintf('DOCLUSTERING - unsupported algorithm "%s"',rep_method))
    end %of switch/case
end
end


function medoid=medoid_SAX(c,clusterNum,SAX_nor_traj,dis_method,varargin)
t=find(c(:,1)==clusterNum);
if isempty(t)
    mmean=[];
else
    dis=Mtx_Distance(SAX_nor_traj(t),SAX_nor_traj(t),'same',dis_method,varargin{:});
    % dis=dis^2;
    sum_dis=sum(dis);
    [s,m]=min(sum_dis);
    mmean=SAX_nor_traj{t(m)};
end
medoid=mmean;
end


function medoid=medoid_SAX_weight(c,clusterNum,SAX_nor_traj,weight,varargin)
t=find(c(:,1)==clusterNum);

if isempty(t)
    mmean=[];
elseif length(t)<3
    mmean=SAX_nor_traj{t(1)};
else
    weight=weight(t);
    for i=1:length(t)
        for j=1:length(t)
            w(i,j)=weight(i)+weight(j);
        end
    end
    dis=Mtx_Distance(SAX_nor_traj(t),SAX_nor_traj(t),'same','Norm',varargin{:});
    dis=dis./w;
    [s,m]=min(sum(dis));
    mmean=SAX_nor_traj{t(m)};
end
medoid=mmean;
end



% function medoid =centre_mediod(c,clusterNum,nor_traj,varargin)
%
% t=find(c(:,1)==clusterNum);
% if isempty (t)
%     medoid=[];
% elseif length(t)<3
%     medoid=nor_traj{t(1)};
% else
%
%     %find distance of objects in cluster
%     dis=Mtx_Distance(nor_traj(t),nor_traj(t),'same','Norm',varargin{:});
%
%     %find the SSE
%     dis=dis.^2;
%     Error=sum(dis);
%     [s,m]=min(Error);
%     medoid=nor_traj{t(m)};
% end
% end

function [medoid inx]=centre_mediod(c,clusterNum,nor_traj_raw,dis,varargin)
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
    inx=t(1);
else
    nor_traj=represent_TS(nor_traj_raw(t),options.rep,varargin{:});
    %find distance of objects in cluster
    
    %dis=Mtx_Distance(nor_traj,nor_traj,'same',varargin{:});
    
    %find the SSE
    dis=dis.^2;
    Error=sum(dis);
    [s,m]=min(Error);
    % medoid=nor_traj{m}; %transfered
    medoid=nor_traj_raw{t(m)};
    inx=t(m);
end
end

function [medoid inx]=centre_mediod_weighted(c,clusterNum,nor_traj_raw,dis,varargin)
%clusterNum
options = struct('rep','RAW','weight',[]);
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
    inx=t(1);
else
    nor_traj=represent_TS(nor_traj_raw(t),options.rep,varargin{:});
    %find distance of objects in cluster
    
    %dis=Mtx_Distance(nor_traj,nor_traj,'same',varargin{:});

    weight=options.weight(t);
    for i=1:length(t)
        for j=1:length(t)
            w(i,j)=weight(i)+weight(j);
        end
    end
    dis=dis./w;

    %find the SSE
    dis=dis.^2;
    Error=sum(dis,2);
    [s,m]=min(Error);
    % medoid=nor_traj{m}; %transfered
    medoid=nor_traj_raw{t(m)};
    inx=t(m);
end
end
