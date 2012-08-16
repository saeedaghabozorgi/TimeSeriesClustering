function [details]=clustering_Hybrid_3Level_anytime3(nor_traj, kk,p,varargin)
plot_show=0;
%% ------------ Level 1-- k-mode --------
disp(['-- START ---------------------------------']);
k=round(length(nor_traj)/30);
if k<kk
    k=kk;
end
parameter1={'l1_dis_method','Euclid','l1_dtw_bound',1,'l1_rep','RAW','l1_alphabet_size',8,'l1_compression_ratio',8,'l1_alg','hier_avg'};
[c1 ,D]=clustering_l1_preclustering(k,p,nor_traj,1,parameter1{:});
aux_details=[2];
details=[1];
i=1;
A=zeros(length(nor_traj),length(nor_traj));
while i<55;
    aux_details=details;
    %% -------------Level 2--CAST--------------------------------------------
    parameter2={'l2_dis_method','DTW','l2_dtw_bound',1,'l2_rep','RAW','l2_alphabet_size',8,'l2_compression_ratio',2};
    [c2 ,error_rate,N_reduction_2lev,D,A]=clustering_l2_purify3(c1,p,nor_traj,D,A,1,parameter2{:});
    [center weight cen_inx]=clustering_l2_making_prototype(c2,p,nor_traj,D,parameter2{:});
    %% --------------Level 3-------------------------------------------------
    parameter3={'l3_dis_method','DTW','l3_dtw_bound',1,'l3_rep','RAW','l3_alphabet_size',8,'l3_compression_ratio',2,'l3_alg','k-medoids-keogh'};
    [c1 ,details,D,A]=clustering_l3_merge(c2,p,k,center,cen_inx',nor_traj,weight,D,A,1,dist_mtx_DTW,parameter3{:});
    calculations_ratio(i)=sum(sum(A))./(length(A)*length(A));
    disp(['  --> D:',num2str(calculations_ratio(i))]);
    %  parameter3={'l3_dis_method','DTW','l3_dtw_bound',1,'l3_rep','RAW','l3_alphabet_size',8,'l3_compression_ratio',2,'l3_alg','k-medoids_weighted'};
    [c3 ,w,v,v]=clustering_l3_merge(c2,p,kk,center,cen_inx',nor_traj,weight,D,A,0,parameter3{:});
    ww(i,:)=w;
    i=i+1;
    disp(['  --',num2str(i),'---------------------------------']);
end
disp(['  --> quality Ultimate:']);
parameter3={'l3_dis_method','DTW','l3_dtw_bound',1,'l3_rep','RAW','l3_alphabet_size',8,'l3_compression_ratio',2,'l3_alg','k-medoids-keogh'};
[c3 ,details]=clustering_l3_merge(c2,p,kk,center,cen_inx',nor_traj,weight,D,A,parameter3{:});

end


