function [c details]=clustering_l3_merge(c,p,k,center,nor_traj,weight,D,A,varargin)
plot_show=0;
options = struct('l3_dis_method','-','l3_dtw_bound',-1,'l3_rep','-','l3_alphabet_size',0,'l3_compression_ratio',0,'l3_alg','-');
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

% to make the prototypes
%return;
disp('level 3');
disp(['  ','DS',':',num2str(length(center)),' | ','K',':',num2str(k),' | ','dis_method',':',options.l3_dis_method,' | ','dtw_bound',':',num2str(options.l3_dtw_bound),' | ','rep',':',options.l3_rep,' | ','alphabet_size',':',num2str(options.l3_alphabet_size),' ','compression_ratio',':',num2str(options.l3_compression_ratio)]);
disp(['  clustering:',num2str(options.l3_alg)]);
l2_clusterCount=max(c);
if l2_clusterCount>k
    if strmatch(options.l3_alg,'hier_avg')
        [c3,Z]=do_Hierarchical_time(center,k,'average',-1,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
    elseif strmatch(options.l3_alg,'k-means')
        [c3,itr]= do_kMeans_time (center,k,'DTW',0,'RAW','dtw_bound',1);
    elseif strmatch(options.l3_alg,'wk-medoids')
        [c3,~]= do_kMediod_time (center,k,0,'weight',weight, 'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
    elseif strmatch(options.l3_alg,'k-medoids')
        [c3,~]= do_kMediod_time (center,k,0,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
    elseif strmatch(options.l3_alg,'k-medoids-keogh')
        nor_traj=represent_TS(center,options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
        Dist=Mtx_Distance(nor_traj,nor_traj,'same','Norm', 'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
        Dist=squareform(Dist);
        [c3,~]= do_kMedoids_keogh(k,Dist); 
    elseif strmatch(options.l3_alg,'k-medoids-keogh-anytime')
        [c3,~]= do_kMediod_time_anytime (center,k,0,c,p,'weight',weight,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
    else
        disp(['  --> Err: please determine the name of algorithm.']);
    end
    for j=1:k
        l3_mems=find(c3==j);
        for i=1:length(l3_mems)
            sub_members=find(c(:,1)==l3_mems(i));
            c(sub_members,2)=j;
        end
    end
else
    c(:,2)=c(:,1);
end
c_center=c(:,1);
c=c(:,2);

%--evaluation-------------------------------------------------------------
if plot_show
    Plot_time_series_luminate(0,0,c3,p,[],center,[],k,0,0.5,4);
    Plot_time_series_luminate(0,0,c,p,[],nor_traj,[],k,0,0.5,5);
end;
[SSEP,SSEC,RI,ARI,purity,BCubed,ConEntropy,f_measure,jacard,FM,NMI,quality]= do_Evaluate(p,c,nor_traj,[],[]);
details=[SSEP,SSEC,RI,ARI,purity,BCubed,ConEntropy,f_measure,jacard,FM,NMI,quality];
disp(['  --> quality:',num2str(quality)]);
end