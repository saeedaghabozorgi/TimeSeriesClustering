function [c error_rate N_reduction_2lev D A]=clustering_l2_purify3(c,p,nor_traj_raw,D,A,plot_show,varargin)
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
%A=zeros(length(nor_traj_raw),length(nor_traj_raw));

for i=1:clusterCount
    newData=find(c(:,1)==i);
    dismatrix=D(newData,newData);
    dis=dismatrix;
    SSE=sum(dis,2)./(size(dis,1)-1);
    info(1,i)=  mean(SSE); %pre_avg_dis
    info(2,i) = std(SSE); %pre_sigma_dis
    info(3,i)= sum(sum(A(newData,newData)))/((length(newData))^2); %chance A

    Nor = dismatrix - min( dismatrix(:) );
    if max( Nor(:) ) ~= 0
        dismatrix = Nor / max( Nor(:) );
    else
        dismatrix=Nor;
    end
    dis=dismatrix;
    sim=1-dis;
    sim(1:length(sim)+1:length(sim)*length(sim))=0;
    Affin=sum(sim,2)./(size(sim,1)-1);
    info(4,i)=  mean(Affin); %pre_avg_sim
    info(5,i) = std(Affin); %pre_sigma_sim
end

[vel,inx]=max(info(1,:).*info(2,:).*(1-info(3,:)));
if sum(1-info(3,:))<0.002
    [vel,inx]=max(info(1,:).*info(2,:));
end

for i=1:clusterCount
    temp_c=[];
    newData=find(c(:,1)==i);
    % representation
    nor_traj=represent_TS(nor_traj_raw,options.l2_rep,varargin{:});
    
    %--------------------------------
   % dist1=Mtx_Distance(nor_traj(newData),nor_traj(newData),'same','org','dis_method',options.l2_dis_method,'dtw_bound',options.l2_dtw_bound,'alphabet_size',options.l2_alphabet_size,'compression_ratio',options.l2_compression_ratio);
    %D(newData,newData)=dist1;
    
    if i~=inx
        temp_c=ones(length(newData),1);
        mu=0;
        sigma=0;
    else
        %--------------------------
      %  -- to check A to not calculate the distance again
%                 for ii=1:length(newData)-1
%                     for jj=ii+1:length(newData)
%                         if A(newData(ii),newData(jj))==0
%                             D(newData(ii),newData(jj))=dis_dtw3(nor_traj{newData(ii)},nor_traj{newData(jj)},length(nor_traj{newData(ii)}));
%                             D(newData(jj),newData(ii))= D(newData(ii),newData(jj));
%                         end
%                     end
%                 end
        %-- to speead up
        dismat=dismat+dismat';
        D(newData,newData)=dismat(newData,newData);
        %----------------------
        
        dist1=D(newData,newData);
        A(newData,newData)=1 ;
        [temp_c,mu,sigma]= do_CAST_time (nor_traj(newData),dist1,-1,'dis_method',options.l2_dis_method,'rep',options.l2_rep,'alphabet_size',options.l2_alphabet_size,'compression_ratio',options.l2_compression_ratio,'dtw_bound',options.l2_dtw_bound);
    end
    c(newData,2)=temp_c;
   % disp(['  --> Pre-cluster#',num2str(i),'  Mems:',num2str(length(newData)),'  Clus:',num2str(max(temp_c)),' avg_sim_l1:(',num2str(info(1,i)),'-',num2str(info(2,i)),')  avg_sim_DTW:',num2str(mu),'-',num2str(sigma),')']);
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