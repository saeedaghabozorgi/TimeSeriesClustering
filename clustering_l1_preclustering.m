function [c distance details]=clustering_l1_preclustering(k,p,nor_traj_raw,plot_show,dist_mtx,varargin)
options = struct('l1_dis_method','-','l1_dtw_bound',0,'l1_rep','SAX','l1_alphabet_size',0,'l1_compression_ratio',0,'l1_alg','-');
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
%% ----------Level 1----------------------------------------------

disp('Stage 1');
disp(['  ','Algo',':',options.l1_alg,' | ','DS',':',num2str(length(nor_traj_raw)),' | ','k',':',num2str(k),' | ','dis_method',':',options.l1_dis_method,' | ','bound',':',num2str(options.l1_dtw_bound),' | ','rep',':',options.l1_rep,' | ','alphabet_size',':',num2str(options.l1_alphabet_size),' | ','compression_ratio',':',num2str(options.l1_compression_ratio)]);
c=[];
%k1=k;

nor_traj=represent_TS(nor_traj_raw,options.l1_rep,'alphabet_size',options.l1_alphabet_size,'compression_ratio',options.l1_compression_ratio);
if ~isempty (dist_mtx)
    distance=dist_mtx_DTW(newData,newData);
else
    distance=Mtx_Distance(nor_traj,nor_traj,'same','Org', 'dis_method',options.l1_dis_method,'dtw_bound',options.l1_dtw_bound,'rep',options.l1_rep,'alphabet_size',options.l1_alphabet_size,'compression_ratio',options.l1_compression_ratio);
end
if strmatch(options.l1_alg,'k-modes')
    [c,itr]= do_kModes_time(nor_traj_raw,k,0,'dis_method',options.l1_dis_method,'dtw_bound',options.l1_dtw_bound,'rep',options.l1_rep,'alphabet_size',options.l1_alphabet_size,'compression_ratio',options.l1_compression_ratio);
elseif strmatch(options.l1_alg,'k-medoids')
    [c,~]= do_kMediod_time (nor_traj_raw,k,0,distance,'dis_method',options.l1_dis_method,'dtw_bound',options.l1_dtw_bound,'rep',options.l1_rep,'alphabet_size',options.l1_alphabet_size,'compression_ratio',options.l1_compression_ratio);
elseif strmatch(options.l1_alg,'hier_avg')
    [c,~,distance]= do_Hierarchical_time (nor_traj_raw,k,'average',-1,distance,'dis_method',options.l1_dis_method,'dtw_bound',options.l1_dtw_bound,'rep',options.l1_rep,'alphabet_size',options.l1_alphabet_size,'compression_ratio',options.l1_compression_ratio);
elseif strmatch(options.l1_alg,'hier_single')
    [c,~,distance]= do_Hierarchical_time (nor_traj_raw,k,'single',-1,distance,'dis_method',options.l1_dis_method,'dtw_bound',options.l1_dtw_bound,'rep',options.l1_rep,'alphabet_size',options.l1_alphabet_size,'compression_ratio',options.l1_compression_ratio);
elseif strmatch(options.l1_alg,'k-medoids-keogh')
    dis=squareform(distance);
    [c,~]= do_kMedoids_keogh(k,dis);
end
[SSEP,SSEC,RI,ARI,purity,ConEntropy,f_measure,jacard,FM,NMI,CSM,quality]= do_Evaluate(p,c,nor_traj_raw,[],[]);
details=[SSEP,SSEC,RI,ARI,purity,ConEntropy,f_measure,jacard,FM,NMI,CSM,quality];
error_rate=Calculate_error_rate(c,p);
disp(['  --> Number of clusters:',num2str(k),' | ','quality:',num2str(quality),' | ','error_rate:', num2str(error_rate)]);

if plot_show
    Plot_time_series_luminate(0,0,c,p,[],nor_traj_raw,[],k,0,0.5,1);
end
end
