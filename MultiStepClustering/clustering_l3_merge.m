function [c details D A]=clustering_l3_merge(c_inp,p,k,center,cen_inx,nor_traj_raw,weight,D,A,plot_show,evaluation,dist_mtx_DTW,varargin)
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
disp('Stage 3');
disp(['  ','DS',':',num2str(length(center)),' | ','K',':',num2str(k),' | ','dis_method',':',options.l3_dis_method,' | ','dtw_bound',':',num2str(options.l3_dtw_bound),' | ','rep',':',options.l3_rep,' | ','alphabet_size',':',num2str(options.l3_alphabet_size),' ','compression_ratio',':',num2str(options.l3_compression_ratio)]);
disp(['  clustering:',num2str(options.l3_alg)]);
l2_clusterCount=max(c_inp);


nor_traj=represent_TS(nor_traj_raw,options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);

if ~isempty (dist_mtx_DTW)
    D(cen_inx,cen_inx)=dist_mtx_DTW(cen_inx,cen_inx);
else
    % to do: check A to not calculate the distance again
    %     for ii=1:length(cen_inx)-1
    %         for jj=ii+1:length(cen_inx)
    %             if A(cen_inx(ii),cen_inx(jj))==0
    %                 D(cen_inx(ii),cen_inx(jj))=dis_dtw3(nor_traj{cen_inx(ii)},nor_traj{cen_inx(jj)},length(nor_traj{cen_inx(ii)}));
    %                 D(cen_inx(jj),cen_inx(ii))= D(cen_inx(ii),cen_inx(jj));
    %             end
    %         end
    %     end
    D(cen_inx,cen_inx)= Mtx_Distance(nor_traj(cen_inx),nor_traj(cen_inx),'same','org','dis_method',options.l3_dis_method,'dtw_bound',options.l3_dtw_bound,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio);
end
A(cen_inx,cen_inx)=1;
Dist=D(cen_inx,cen_inx);
if l2_clusterCount<k
    disp(['  --> Err: stage2 clusters is less than k.']);
    return;
end

if strmatch(options.l3_alg,'hier_avg')
    [c3,Z]=do_Hierarchical_time(center,k,'average',-1,Dist,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
elseif strmatch(options.l3_alg,'hier_single')
    [c3,Z]=do_Hierarchical_time(center,k,'single',-1,Dist,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
elseif strmatch(options.l3_alg,'k-means')
    [c3,itr]= do_kMeans_time (center,k,'DTW',0,'RAW','dtw_bound',1);
elseif strmatch(options.l3_alg,'k-medoids_weighted')
    [c3,~]= do_kMediod_time (center,k,0,Dist,'weight',weight, 'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
elseif strmatch(options.l3_alg,'k-medoids')
    [c3,~]= do_kMediod_time (center,k,0,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
elseif strmatch(options.l3_alg,'k-medoids-keogh')
    Dist=squareform(Dist);
    [c3,~]= do_kMedoids_keogh(k,Dist);
elseif strmatch(options.l3_alg,'k-medoids-keogh-anytime')
    [c3,~]= do_kMediod_time_anytime (center,k,0,c,p,'weight',weight,'dis_method',options.l3_dis_method,'rep',options.l3_rep,'alphabet_size',options.l3_alphabet_size,'compression_ratio',options.l3_compression_ratio,'dtw_bound',options.l3_dtw_bound);
else
    disp(['  --> Err: please determine the name of algorithm.']);
end



c=labeling(c3,c_inp,k);

%--evaluation-------------------------------------------------------------

if plot_show
    % to plot the prototypes
    %  Plot_time_series_luminate(0,0,c3,p(cen_inx),[],center,[],k,2,0.5,4);
    %     Plot_time_series_luminate(0,0,c3,p,[],center,[],k,0,0.5,4);
    Plot_time_series_luminate(0,0,c,p,[],nor_traj_raw,[],k,0,0.5,5);
end;
details=[];
if evaluation
    [SSEP,SSEC,RI,ARI,purity,ConEntropy,f_measure,jacard,FM,NMI,CSM,quality]= do_Evaluate(p,c,nor_traj_raw,[],[]);
    details=[SSEP,SSEC,RI,ARI,purity,ConEntropy,f_measure,jacard,FM,NMI,CSM,quality];
    disp(['  --> quality:',num2str(quality)]);
end
end

function c=labeling(c3,c_inp,k)
for j=1:k
    l3_mems=find(c3==j);
    for i=1:length(l3_mems)
        sub_members=find(c_inp(:,1)==l3_mems(i));
        c_inp(sub_members,2)=j;
    end
end
c=c_inp(:,2);
end