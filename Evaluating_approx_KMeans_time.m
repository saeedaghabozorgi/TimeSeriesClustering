function Evaluating_approx_KMeans_time

fnames = dir('..\data\dataset UCR\All train\*');
tradata = cell(1,length(fnames));

for k=3:length(fnames)
    
    fname = fnames(k).name;
    files_name{k-2}=fname;
end

for dataset_no=7:  length(files_name)
    file_name=['..\data\dataset UCR\All train\' files_name{dataset_no}];
    disp(file_name);
    train_data = importdata(file_name);
    TRAIN_class_labels = train_data(:,1);     % Pull out the class labels.
    k=length(unique(TRAIN_class_labels));
    rows=size(train_data,1);
    data_len= size(train_data,2)-1;
    compression_ratio=8;
    alphabet_size = 6;
    
    data_len=(floor(data_len/compression_ratio))*compression_ratio;
    [nor_traj,t_traj]=Import_Data_UCR(1,rows,2,data_len+1,file_name);
    
    nseg = data_len/compression_ratio;
    SAX_nor_traj=[];
    for i=1:length(nor_traj)
        SAX_nor_traj{i}= timeseries2symbol(nor_traj{i}, data_len, nseg, alphabet_size);
    end
    
    %----------Level 1----------------------------------------------
    all_details=[];
    for test_no=1:3
        for k1=2:k
            k1
            c=[];
            %----------------------------------------------------------
            
            temp=ones(rows,1);   % initialize as zero vector
            center_SAX=[];
            p = randperm(rows);      % random initialization
            for i=1:k1
                center_SAX{i}=SAX_nor_traj{p(i)};
                temp(p(i))=i;
            end
            itr=1;
            while 1,
                dis=Mtx_SAX_Max_Distance(SAX_nor_traj,center_SAX,'',alphabet_size,compression_ratio);
                
                [z,c]=min(dis,[],2);  % find group matrix g
                
                if (c==temp | itr==100),
                    break;          % stop the iteration
                else
                    temp=c;         % copy group matrix to temporary variable
                end
                for i=1:k1
                    center{i}=centre_mean(c,i,nor_traj);
                    center_SAX{i}= timeseries2symbol(center{i}, data_len, nseg, alphabet_size);
                end
                itr=itr+1;
            end
            
            clusterCount=max(c);
            %----------------------------------------------------------
            
            cp=[];
            [cp,type]= do_kMeans_time (nor_traj,k1,'Euclid',0)  ;
            %----------------------------------------------------------

            
            RS= calculate_Cluster_Rand_statistic(c,cp);
            purity= Calculate_Cluster_Purity(c,cp,2);
            [BCubed,f_precision,f_recall]=Calculate_Cluster_BCubed_precision(c,cp);
            entropy=Calculate_cluster_entropy(c,cp);
            fm=Calculate_Cluster_F_measure(c,cp);
            
            details(k1,1)=k;
            details(k1,2)=rows;
            details(k1,3)=data_len;
            details(k1,6)=RS;
            details(k1,7)=purity;
            details(k1,8)=BCubed;
            details(k1,9)=f_precision;
            details(k1,10)=f_recall;
            details(k1,11)=entropy;
            details(k1,12)=fm;
            
            % Import the file
            %     p=train_data(:,1);
            %     if(min(p)==0)
            %         p=p+1;
            %     end
            %
            %     RS= calculate_Cluster_Rand_statistic(c,p);
            %     purity= Calculate_Cluster_Purity(c,p,2);
            %     sim= Calculate_Cluster_Similarity(c,p);
        end
        if isempty(all_details)
            all_details=details;
        else
            
            all_details=all_details+details;
        end
    end
    all_details=all_details/test_no;
end
end

function mmean=centre_mean(c,clusterNum,nor_traj)
t=find(c(:,1)==clusterNum);
if isempty(t)
    p = randperm(length(nor_traj));      % random initialization
    mmean= nor_traj{p(1)};
else
    n=length(nor_traj{1});
    cluster_mem=zeros(1,n);
    for j=1:size(t,1)
        cluster_mem(j,:)=nor_traj{t(j)};
    end
    mmean=mean(cluster_mem,1);
end
end