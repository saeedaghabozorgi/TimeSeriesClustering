function c=clustering_l1_preclustering(k,p,nor_traj,varargin)
plot_show=0;
options = struct('l1_dis_method','SAXminDis','l1_dtw_bound',1,'l1_rep','SAX','l1_alphabet_size',8,'l1_compression_ratio',8);
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
k1=round(length(nor_traj)/30);
if k1<k
    k1=k;
end
disp('level 1');
disp(['  ','DS',':',num2str(length(nor_traj)),' | ','k',':',num2str(k1),' | ','dis_method',':',options.l1_dis_method,' | ','rep',':',options.l1_rep,' | ','alphabet_size',':',num2str(options.l1_alphabet_size),' | ','compression_ratio',':',num2str(options.l1_compression_ratio)]);
c=[];
%k1=k;
[c,itr]= do_kModes_time(nor_traj,k1,0,'dis_method',options.l1_dis_method,'rep',options.l1_rep,'alphabet_size',options.l1_alphabet_size,'compression_ratio',options.l1_compression_ratio);
disp(['  --> DS:',num2str(length(nor_traj)),'  clusters:',num2str(k1)]);
%   details_l1=[k,SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality];
if plot_show
    Plot_time_series_luminate(0,0,c,p,[],nor_traj,[],k1,0,0.5,1);
    [SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality]= do_Evaluate(p,c,nor_traj,[],[]);
    disp(['  --> quality:',num2str(quality)]);
end
end
