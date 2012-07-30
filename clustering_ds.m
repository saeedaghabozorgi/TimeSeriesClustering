function clustering_ds
file_name=['..\Datasets\Face(all)\FaceAll.mat'];
train_data = importdata(file_name);
lable_file_name=['..\Datasets\Face(all)\Label_FaceAll.mat'];
train_data_lable = importdata(lable_file_name);
%-----------to min
d_data=[];
lab=[];
for i=1:14
    x=find(train_data_lable==i);
    x=x(1:5);
    d_data=[d_data;train_data(x,:)];
    lab=[lab;train_data_lable(x,:)];
end
train_data=d_data;
train_data_lable=lab;
%-----------------

% train_data_lable=train_data_lable(p,:);
p=train_data_lable(:,1);
k=length(unique(p));
rows=size(train_data,1);
data_len= size(train_data,2);
[nor_traj,~]=NormalizeTS(train_data);

%% clustering
%-----------------------------
details(1,:)=clustering_Hybrid_3Level(nor_traj,k,p);
%-----------------------------
% parameter1={'l1_dis_method','DTW','l1_dtw_bound',1,'l1_rep','RAW','l1_alphabet_size',8,'l1_compression_ratio',2,'l1_alg','k-medoids-keogh'};
% c=clustering_l1_preclustering(k,p,nor_traj,parameter1{:});
end
