function clustering_ds
%file_name=['..\dataset UCR large\Face(all)\FaceAll.mat'];
%file_name=('..\dataset CBF\ds_CBF_300.mat');
file_name=('..\dataset CC\ds_CC_300.mat');
train_data = importdata(file_name);
%lable_file_name=['..\dataset UCR large\prepared_k-medoids_label\ds_FaceAll_label.mat'];
%lable_file_name=['..\dataset CBF\ds_CBF_300_label.mat'];
lable_file_name=['..\dataset CC\ds_CBF_300_kmedoids_DTW_label.mat'];
train_data_lable = importdata(lable_file_name);
%-----------to min
% d_data=[];
% lab=[];
% for i=1:14
%     x=find(train_data_lable==i);
%     x=x(1:2);
%     d_data=[d_data;train_data(x,:)];
%     lab=[lab;train_data_lable(x,:)];
% end
% train_data=d_data;
% train_data_lable=lab;
%-----------------

% train_data_lable=train_data_lable(p,:);
p=train_data_lable(:,1);
k=length(unique(p));
rows=size(train_data,1);
data_len= size(train_data,2);
[nor_traj,~]=NormalizeTS(train_data);

%% clustering
%-----------------------------
%save_k_medoids_lables(k);
%paralele_DTW(nor_traj);
%-----------------------------
 load('..\dataset CC\prepared_mtx_CC\dismat_DTW_CC_300.mat','dismat');
 dismat=dismat+dismat';
% dismat=squareform(dismat);
details(1,:)=clustering_Hybrid_3Level(nor_traj,k,p,dismat);
%-----------------------------
%details(1,:)=clustering_Hybrid_3Level_anytime3(nor_traj,k,p);
%-----------------------------
%  parameter1={'l1_dis_method','DTW','l1_dtw_bound',1,'l1_rep','RAW','l1_alphabet_size',8,'l1_compression_ratio',2,'l1_alg','k-medoids-keogh'};
%  c=clustering_l1_preclustering(k,p,nor_traj,parameter1{:});
%-----------------------------
%  load('..\dataset CBF\prepared_mtx_CBF\dismat_DTW_CBF_300.mat','dismat');
%  dismat=dismat+dismat';
%  dismat=squareform(dismat);
% [c,~]= do_kMedoids_keogh(k,dismat);
% [SSEP,SSEC,RI,ARI,purity,BCubed,ConEntropy,f_measure,jacard,FM,NMI,quality]= do_Evaluate(p,c,nor_traj,[],[]);
% details=[SSEP,SSEC,RI,ARI,purity,BCubed,ConEntropy,f_measure,jacard,FM,NMI,quality];
%disp(['  --> quality:',num2str(quality)]);
end

function save_k_medoids_lables(k)
        load('..\dataset UCR large\prepared_mtx_UCR_large\dismat_DTW_FaceAll.mat','dismat');
        D=dismat+dismat';
        D=squareform(D);
        [label,~]= do_kMedoids_keogh(k,D);
        save('..\dataset UCR large\prepared_k-medoids_label\ds_FaceAll_label.mat', 'label')
end

function paralele_DTW(nor_traj)
matlabpool open 2
poolSize = matlabpool('size');
if poolSize == 0
   error('parallel:demo:poolClosed', ...
        'This demo needs an open MATLAB pool to run.');
end
fprintf('This demo is running on %d MATLABPOOL workers.\n', ...
    matlabpool('size'));

x=nor_traj;
data_n = length(x);
dismat=zeros(data_n,data_n);
parfor  i = 1:data_n
    disp(num2str(i));
    dismat(i,:)=rowcalc(i,data_n,x);
end
save('dismat_DTW_FaceAll.mat', 'dismat')
 matlabpool close
end

function dismatrix=rowcalc(row,data_n,x)
dismatrix=zeros(1,data_n);
    a=x{row};
    for j = row:data_n,
        if row ~= j
            b=x{j};
            dis=dis_dtw3(a,b,size(a,2));
            dismatrix(1, j) =dis;
        end
    end
end
