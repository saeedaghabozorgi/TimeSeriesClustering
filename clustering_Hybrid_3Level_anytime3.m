function ww=clustering_Hybrid_3Level_anytime3(nor_traj, kk,p,dist_mtx_DTW,varargin)
plot_show=0;
%% ------------ Level 1-- k-mode --------
disp(['-- START ---------------------------------']);
k=round(length(nor_traj)/50);
if k<kk
    k=kk;
end
parameter1={'l1_dis_method','Euclid','l1_dtw_bound',1,'l1_rep','RAW','l1_alphabet_size',8,'l1_compression_ratio',4,'l1_alg','k-medoids-keogh'};
[c1 ,D]=clustering_l1_preclustering(k,p,nor_traj,0,[],parameter1{:});
aux_details=[2];
details=[1];
i=1;
A=zeros(length(nor_traj),length(nor_traj));
while i<50;
    aux_details=details;
    %% -------------Level 2--CAST--------------------------------------------
    parameter2={'l2_dis_method','DTW','l2_dtw_bound',1,'l2_rep','RAW','l2_alphabet_size',8,'l2_compression_ratio',2};
    [c2 ,error_rate,N_reduction_2lev,D,A]=clustering_l2_purify(c1,p,nor_traj,D,A,0,dist_mtx_DTW,parameter2{:});
    [center weight cen_inx]=clustering_l2_making_prototype(c2,p,nor_traj,D,0,parameter2{:});
    %% --------------Level 3-------------------------------------------------
    parameter3={'l3_dis_method','DTW','l3_dtw_bound',1,'l3_rep','RAW','l3_alphabet_size',8,'l3_compression_ratio',2,'l3_alg','k-medoids-keogh'};
    [c1 ,details,D,A]=clustering_l3_merge(c2,p,k,center,cen_inx',nor_traj,weight,D,A,0,dist_mtx_DTW,parameter3{:});
    calculations_ratio(i)=sum(sum(A))./(length(A)*length(A));
    disp(['  --> D:',num2str(calculations_ratio(i))]);
     parameter3={'l3_dis_method','DTW','l3_dtw_bound',1,'l3_rep','RAW','l3_alphabet_size',8,'l3_compression_ratio',2,'l3_alg','k-medoids-keogh'};
    [~ ,w,~,~]=clustering_l3_merge(c2,p,kk,center,cen_inx',nor_traj,weight,D,A,0,dist_mtx_DTW,parameter3{:});
    
%     [ct,~]= do_kMedoids_keogh(k,squareform(D));
%     [SSEP,SSEC,RI,ARI,purity,ConEntropy,f_measure,jacard,FM,NMI,CSM,quality]= do_Evaluate(p,ct,nor_traj,[],[]);
%     det(i,:)=[SSEP,SSEC,RI,ARI,purity,ConEntropy,f_measure,jacard,FM,NMI,CSM,quality];
    
    ww(i,:)=w;
    i=i+1;
    disp(['  --',num2str(i),'---------------------------------']);
end
end


