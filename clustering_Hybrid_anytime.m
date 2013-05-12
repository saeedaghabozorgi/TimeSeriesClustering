function clustering_Hybrid_anytime()

for dataset_no=1:1
    train_data = CBF_Generator(10*dataset_no);
    TRAIN_class_labels = train_data(:,1);     % Pull out the class labels.
    p=train_data(:,1);
    k=length(unique(TRAIN_class_labels));
    rows=size(train_data,1);
    data_len= size(train_data,2)-1;
    [nor_traj,~]=Import_Data_UCR1(1,rows,2,data_len+1,train_data);
    cluster_count(dataset_no)=k;
    pp{dataset_no}=p;
    ds{dataset_no}=nor_traj;
    dataset_no
    %parameter={'l1_dis_method','Euclid','l1_dtw_bound',0.5,'l1_rep','SAX','l1_alphabet_size',8,'l1_compression_ratio',6};
    %   parameter=[parameter,  'l2_dis_method','Euclid','l2_dtw_bound',.2,'l2_rep','RAW','l2_alphabet_size',8,'l2_compression_ratio',1];
    %parameter=[parameter,  'l3_dis_method','DTW','l3_dtw_bound',0.3,'l3_rep','RAW','l3_alphabet_size',8,'l3_compression_ratio',6,'l3_alg','k-medoid'];
    %   details(dataset_no,:)=clustering_Hybrid_3Level_anytime(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no},parameter{:});
    nor_traj_raw=nor_traj;
    for itr=1 : 9
        itr
        nor_traj=represent_TS(nor_traj_raw,'SAX','alphabet_size',8,'compression_ratio',10-itr);
        dismatrix=Mtx_Distance(nor_traj,nor_traj,'same','dis_method','DTW','dtw_bound',0.9);
        dismatrix=squareform(dismatrix);
        [c, cost, medoids] = do_kMedoids_keogh(k,dismatrix)
        % [c,itr]= do_kModes_time(nor_traj,k1,0,parameter{:});
        
        [SSEP,SSEC,RI,ARI,purity,BCubed,ConEntropy,f_measure,jacard,FM,NMI,quality]= do_Evaluate(p,c,nor_traj,[],[]);
        details(itr,:)=[SSEP,SSEC,RI,ARI,purity,BCubed,ConEntropy,f_measure,jacard,FM,NMI,quality];
    end
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