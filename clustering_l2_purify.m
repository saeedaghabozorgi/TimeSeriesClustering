function [c error_rate N_reduction_2lev D A]=clustering_l2_purify(c,p,nor_traj_raw,D,A,plot_show,dist_mtx_DTW,varargin)

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
disp('Stage 2');
disp(['  ','DS',':',num2str(length(c)),' | ','dis_method',':',options.l2_dis_method,' | ','dtw_bound',':',num2str(options.l2_dtw_bound),' | ','rep',':',options.l2_rep,' | ','alphabet_size',':',num2str(options.l2_alphabet_size),' ','compression_ratio',':',num2str(options.l2_compression_ratio)]);
clusterCount=max(c);
%A=zeros(length(nor_traj_raw),length(nor_traj_raw));

for i=1:clusterCount
    temp_c=[];
    newData=find(c(:,1)==i);
    
    %-------------------
    % only for print
    dismatrix=D(newData,newData);
    Nor = dismatrix - min( dismatrix(:) );
    if max( Nor(:) ) ~= 0
        dismatrix = Nor / max( Nor(:) );
    else
        dismatrix=Nor;
    end
    dis=dismatrix;
    sim=1-dis;
    sim(1:length(sim)+1:length(sim)*length(sim))=0;
    
    Affin1=sum(sim,2)./(size(sim,1)-1);
    pre_avg_sim=  mean(Affin1);
    pre_sigma_sim = std(Affin1);
    %        figure;
    %        dd=squareform(sim);
    %     hist(dd);
    %--------------
    % representation
    nor_traj=represent_TS(nor_traj_raw,options.l2_rep,'alphabet_size',options.l2_alphabet_size,'compression_ratio',options.l2_compression_ratio);
    %--------------------------
    %--------------------------------
    %dist1=
    if ~isempty (dist_mtx_DTW)
        D(newData,newData)=dist_mtx_DTW(newData,newData);
    else
        D(newData,newData)= Mtx_Distance(nor_traj(newData),nor_traj(newData),'same','org','dis_method',options.l2_dis_method,'dtw_bound',options.l2_dtw_bound,'alphabet_size',options.l2_alphabet_size,'compression_ratio',options.l2_compression_ratio);
     
    end
    A(newData,newData)=1 ;
    %-------------------
    % only for print
    dismatrix=D(newData,newData);
    Nor = dismatrix - min( dismatrix(:) );
    if max( Nor(:) ) ~= 0
        dismatrix = Nor / max( Nor(:) );
    else
        dismatrix=Nor;
    end
    dis=dismatrix;
    sim=1-dis;
    sim(1:length(sim)+1:length(sim)*length(sim))=0;
    Affin2=sum(sim,2)./(size(sim,1)-1);
    DTW_avg_sim=  mean(Affin2);
    DTW_sigma_sim = std(Affin2);
    %     figure;
    %       hist(Affin2);
    %--------------
    
    dist1=D(newData,newData);
    if length(newData)<5
        temp_c=ones(length(newData),1);
        mu=0;
        sigma=0;
    else
        [temp_c]= do_CAST_time (nor_traj(newData),dist1,min(Affin1),'dis_method',options.l2_dis_method,'rep',options.l2_rep,'alphabet_size',options.l2_alphabet_size,'compression_ratio',options.l2_compression_ratio,'dtw_bound',options.l2_dtw_bound);
    end
    c(newData,2)=temp_c;
    disp(['  --> Pre-cluster#',num2str(i),'  Mems:',num2str(length(newData)),'  Clus:',num2str(max(temp_c)),' avg_sim_l1:(',num2str(pre_avg_sim),'-',num2str(pre_sigma_sim),')  avg_sim_DTW:',num2str(DTW_avg_sim),'-',num2str(DTW_sigma_sim),')']);
end
%-------
% to map raw objects to new clusters
c(:,3)=c(:,1)*100000+c(:,2);
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