function [center weight]=clustering_l2_making_prototype(c,nor_traj,varargin)
plot_show=0;
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
l2_clusterCount=max(c(:,4));
for i=1:l2_clusterCount;
    %center{i}=centre_mediod(c(:,4),i,nor_traj,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
    center{i}=centre_randoid(c(:,4),i,nor_traj,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
    weight(i,1)=length(find(c(:,4)==i));
end
if plot_show 
    % Plot_time_series_luminate(0,0,c(:,4),p,center,nor_traj,[],l2_clusterCount,2,0.2,2); 
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
