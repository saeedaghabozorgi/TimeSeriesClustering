function [center weight cen_inx]=clustering_l2_making_prototype(c,p,nor_traj,dist,plot_show,varargin)
options = struct('l3_dis_method','SAXminDis','l3_dtw_bound',1,'l3_rep','SAX','l3_alphabet_size',8,'l3_compression_ratio',8,'l3_alg','k-medoid');
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
disp('  Making prototype ...');
center=[];
l2_clusterCount=max(c);
for i=1:l2_clusterCount;
    newData=find(c(:,1)==i);
    dis=dist(newData,newData);
    [cent inx]=centre_mediod(c,i,nor_traj,dis,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
    center{i}=cent;
    cen_inx(i)=inx;
   % center{i}=centre_randoid(c,i,nor_traj,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
    weight(i,1)=length(find(c==i));
end
if plot_show 
    %Plot_time_series_luminate(0,0,c,p,center,nor_traj,[],l2_clusterCount,2,0.2,3); 
    c=[1:1:l2_clusterCount]';
    p=p(cen_inx);
    Plot_time_series_luminate(0,0,c,p,[],nor_traj(cen_inx),[],l2_clusterCount,0,0.8,3); 
end;
end

function randoid=centre_randoid(c,clusterNum,nor_traj_raw,varargin)
%clusterNum
t=find(c(:,1)==clusterNum);
if isempty (t)
    randoid=[];
else
    randoid=nor_traj_raw{t(1)};
end
end

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
