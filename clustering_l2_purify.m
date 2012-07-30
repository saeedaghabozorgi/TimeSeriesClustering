function [c error_rate N_reduction_2lev D A]=clustering_l2_purify(c,p,nor_traj_raw,D,varargin)
plot_show=0;
options = struct('l2_dis_method','SAXminDis','l2_dtw_bound',1,'l2_rep','SAX','l2_alphabet_size',8,'l2_compression_ratio',8);
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
%% ----------Level 2--CAST--------------------------------------------
disp('level 2');
disp(['  ','DS',':',num2str(length(c)),' | ','dis_method',':',options.l2_dis_method,' | ','dtw_bound',':',num2str(options.l2_dtw_bound),' | ','rep',':',options.l2_rep,' | ','alphabet_size',':',num2str(options.l2_alphabet_size),' ','compression_ratio',':',num2str(options.l2_compression_ratio)]);
clusterCount=max(c);
A=zeros(length(nor_traj_raw),length(nor_traj_raw));
for i=1:clusterCount
    temp_c=[];
    newData=find(c(:,1)==i);
    % representation
    nor_traj=represent_TS(nor_traj_raw,options.l2_rep,varargin{:});
    dist1=Mtx_Distance(nor_traj(newData),nor_traj(newData),'same','org','dis_method',options.l2_dis_method,'dtw_bound',options.l2_dtw_bound,'alphabet_size',options.l2_alphabet_size,'compression_ratio',options.l2_compression_ratio);
    D(newData,newData)=dist1;
    A(newData,newData)=1 ;
    if length(newData)<5
        temp_c=ones(length(newData),1);
    else
        [temp_c]= do_CAST_time (nor_traj(newData),dist1,-1,'dis_method',options.l2_dis_method,'rep',options.l2_rep,'alphabet_size',options.l2_alphabet_size,'compression_ratio',options.l2_compression_ratio,'dtw_bound',options.l2_dtw_bound);
    end
    c(newData,2)=temp_c;
    disp(['  --> cluster:',num2str(i),'  Mems:',num2str(length(newData)),'  Clus:',num2str(max(temp_c))]);
end
%-------
% to map raw objects to new clusters
c(:,3)=c(:,1)*10000+c(:,2);
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

%--evaluation---------------------------------------------------------
if plot_show
    Plot_time_series_luminate(0,0,c(:,4),p,[],nor_traj_raw,[],l2_clusterCount,2,0.5,2);
end;
c=c(:,4);
purity2=Calculate_Cluster_Purity(c,p,1);
qual_2lev=Calculate_Cluster_correct_ratio(c,p);
N_reduction_2lev=1-l2_clusterCount/length(nor_traj);
error_rate=Calculate_error_rate(c,p);
%[SSEP,SSEC,RI,ARI,purity,BCubed,ConEntropy,f_measure,jacard,FM,NMI,quality]= do_Evaluate(p,c,nor_traj,[],[]);
disp(['  --> Number of clusters:',num2str(l2_clusterCount),' | ','error_rate:', num2str(error_rate),' | ', 'reduction:', num2str(N_reduction_2lev),' | ','correct_rate:', num2str(qual_2lev),' | ','purity:', num2str(purity2)]);


end