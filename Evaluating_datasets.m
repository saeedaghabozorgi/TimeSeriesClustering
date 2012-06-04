function Evaluating_datasets

fnames = dir('..\data\dataset UCR\All train\*');
tradata = cell(1,length(fnames));

for k=3:length(fnames)
    
    fname = fnames(k).name;
    files_name{k-2}=fname;
end


for dataset_no=1:length(files_name)
    file_name=['..\data\dataset UCR\All train\' files_name{dataset_no}];
    disp(file_name);
    train_data = importdata(file_name);
    TRAIN_class_labels = train_data(:,1);     % Pull out the class labels.
    p=train_data(:,1);
    if(min(p)==-1)
        p(p==-1)=2;
    end
    if(min(p)==0)
        p(p==0)=2;
    end
      
    p_clust_count=[];
    k=length(unique(TRAIN_class_labels));
   
    rows=size(train_data,1);
    data_len= size(train_data,2)-1;
    
    
    %----------------------------------------
    %   fprintf('Dimensionaly reduction of  data .. \n\r');
    
    compression_ratio=8;
    alphabet_size = 8;
    

    data_len=(floor(data_len/compression_ratio))*compression_ratio;
    [nor_traj,t_traj]=Import_Data_UCR(1,rows,2,data_len+1,file_name);
    
    nseg = data_len/compression_ratio;

    b=unique(p);
    counts = histc(p, b);
    counts=counts/rows;
 
    avg=mean(counts);
    variance=var(counts);
    stdev=std(counts);
    
   
    details(dataset_no,:)=[k,rows,data_len,avg,variance,stdev]
end
end
