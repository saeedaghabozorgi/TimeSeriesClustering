function details=clustering_Hybrid_3Level_anytime(nor_traj, k,p,varargin)
%details(dataset_no,:)=clustering_Hybrid_3Level_anytime(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no},parameter{:});
%% ------------ Level 1-- k-mode --------
parameter1={'l1_dis_method','Euclid','l1_dtw_bound',1,'l1_rep','SAX','l1_alphabet_size',8,'l1_compression_ratio',4};
c1=clustering_l1_preclustering(k,p,nor_traj,parameter1{:});
itr=1;
while itr<8
    disp(['itr:',num2str(itr),' --------------------------------------------------']);
    %% ----------Level 2--CAST--------------------------------------------
    parameter2={'l2_dis_method','DTW','l2_dtw_bound',.1+itr*0.1,'l2_rep','SAX','l2_alphabet_size',8,'l2_compression_ratio',9-itr};
    [c2,error_rate_2lev,N_reduction_2lev]=clustering_l2_purify(c1,p,nor_traj,parameter2{:});
    [center weight]=clustering_l2_making_prototype(c2,nor_traj,parameter2{:});
    %% ----------Level 3-------------------------------------------------
    k3=max(c1);
    parameter3={'l3_dis_method','DTW','l3_dtw_bound',.1+itr*0.1,'l3_rep','SAX','l3_alphabet_size',8,'l3_compression_ratio',9-itr,'l3_alg','k-medoid'};
    [c3,det3]=clustering_l3_merge(c2,p,k3,center,nor_traj,weight,parameter3{:});
    details(itr,:)=[det3,error_rate_2lev];
    c1=c3;
    itr=itr+1;
end
parameter4={'l3_dis_method','DTW','l3_dtw_bound',.9,'l3_rep','RAW','l3_alphabet_size',8,'l3_compression_ratio',9-itr,'l3_alg','k-medoid'};
c=clustering_l3_merge(c2,p,k,center,nor_traj,weight,parameter4{:});
[SSEP,SSEC,RI,ARI,purity,BCubed,ConEntropy,f_measure,jacard,FM,NMI,quality]= do_Evaluate(p,c,nor_traj,[],[]);
details(itr,:)=[SSEP,SSEC,RI,ARI,purity,BCubed,ConEntropy,f_measure,jacard,FM,NMI,quality,error_rate_2lev];
end

