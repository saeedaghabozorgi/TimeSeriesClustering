function clustering_ds_CC_online
detl=[];
for dataset_no=7:10
    train_data = ds_Generator_CC(10*dataset_no);
    TRAIN_class_labels = train_data(:,1);     % Pull out the class labels.
    p=train_data(:,1);
    k=length(unique(TRAIN_class_labels));
    rows=size(train_data,1);
    data_len= size(train_data,2)-1;
    [nor_traj,~]=Import_Data_UCR1(1,rows,2,data_len+1,train_data);
    cluster_count(dataset_no)=k;
    pp{dataset_no}=p;
    ds{dataset_no}=nor_traj;
    %%
    %--- to save CBF and matrix data and lables
%        ds_CBF=train_data(1:300,2:129);
%         foldpath=['..\dataset CC'];
%         file_name=['ds_CC_',int2str(rows)];
%         mkdir([foldpath,'\',file_name]);
%         save([foldpath,'\',file_name,'\',file_name,'.mat'],'train_data');
%         calculate_DTW_matrix_paralele(ds{dataset_no},[foldpath,'\',file_name,'\',file_name,'_dismat_DTW.mat']);
%         load('..\dataset CBF\prepared_mtx_CBF\dismat_DTW_CBF_300.mat','dismat');
%         D=dismat+dismat';
%         D=squareform(D);
%         [label,~]= do_kMedoids_keogh(k,D);
%         save('..\dataset CBF\ds_CBF_300_label.mat', 'label')
    %---------------------------------------
    
        details{dataset_no}=clustering_anytime2(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no},[]);
    
   %    details(dataset_no,:)=clustering_Hybrid_3Level(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no},[]);
    %   details(dataset_no,:)=clustering_Hybrid_3Level_anytime3(ds{dataset_no},cluster_count(dataset_no),pp{dataset_no});
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

function paralele_DTW(nor_traj,filename)
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
save(filename, 'dismat')
matlabpool close
end


function calculate_DTW_matrix_paralele(nor_traj,path)
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
dismat=dismat+dismat';
save(path, 'dismat')
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
