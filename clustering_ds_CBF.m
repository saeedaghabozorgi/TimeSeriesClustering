function clustering_ds_CBF
detl=[];
for dataset_no=10:10
    plot_show=1;
    train_data = CBF_Generator(30*dataset_no);
    TRAIN_class_labels = train_data(:,1);     % Pull out the class labels.
    p=train_data(:,1);
    k=length(unique(TRAIN_class_labels));
    rows=size(train_data,1);
    data_len= size(train_data,2)-1;
    [nor_traj,~]=Import_Data_UCR1(1,rows,2,data_len+1,train_data);
    cluster_count(dataset_no)=k;
    pp{dataset_no}=p;
    ds{dataset_no}=nor_traj;
    dataset_no;
    
     parameter={'l1_dis_method','SAXAPX','l1_dtw_bound',1,'l1_rep','SAX','l1_alphabet_size',8,'l1_compression_ratio',6};
     parameter=[parameter,  'l2_dis_method','DTW','l2_dtw_bound',.07,'l2_rep','RAW','l2_alphabet_size',8,'l2_compression_ratio',6];
     parameter=[parameter,  'l3_dis_method','DTW','l3_dtw_bound',0.7,'l3_rep','SAX','l3_alphabet_size',8,'l3_compression_ratio',6,'l3_alg','k-medoid'];
     details(dataset_no,:)=clustering_Hybrid_3Level(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no},parameter{:});
    
    
    %%    
%     parameter={'l1_dis_method','SAXAPX','l1_dtw_bound',1,'l1_rep','SAX','l1_alphabet_size',8,'l1_compression_ratio',6};
%      parameter=[parameter,  'l2_dis_method','DTW','l2_dtw_bound',.07,'l2_rep','RAW','l2_alphabet_size',8,'l2_compression_ratio',6];
%      parameter=[parameter,  'l3_dis_method','DTW','l3_dtw_bound',0.7,'l3_rep','SAX','l3_alphabet_size',8,'l3_compression_ratio',6,'l3_alg','k-medoid'];
%     details{dataset_no}=clustering_Hybrid_3Level_anytime(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no});
%    detl=[detl;details{dataset_no}];
end
end

function [nor_traj,t_traj]=Import_Data_UCR1(FromUser,ToUser,SinceTime,ToTime,data)
org_traj=[];  % orginal traj
nor_traj=[]; % normalized orginal traj
t_traj=[];

for z=FromUser:ToUser
    
    TS=SinceTime;
    TE=ToTime;
    a=data(z,TS:TE);
    an =(a-mean(a))/std(a);
    inx=length(org_traj)+1;
    nor_traj{inx}=an;
    org_traj{inx}=a;
    t_traj{inx}=(1:1:TE-TS+1);
end
end