function S_plot_sample_class()
    % Saeed, to draw samples of time series from a data set
    clear;
    clc;
    file_name='D:\My Docments\My Research\Research Works\My research_Time series Clustering\Data\dataset UCR\Adiac\Adiac_TEST';
    train_data = importdata(file_name);

    %----labels from data------
    TRAIN_class_labels = train_data(:,1);
    if(min(TRAIN_class_labels)==-1)
        TRAIN_class_labels(TRAIN_class_labels==-1)=max(TRAIN_class_labels)+1;
    end 
    if(min(TRAIN_class_labels)==0)
        TRAIN_class_labels(TRAIN_class_labels==0)=max(TRAIN_class_labels)+1;
    end

    p=TRAIN_class_labels;
    k=length(unique(TRAIN_class_labels));
    rows=size(train_data,1);
    data_len= size(train_data,2)-1;
    [nor_traj_raw,~]=Import_Data_UCR(1,rows,2,data_len+1,file_name);
    
    dist_mtx_file='D:\My Docments\My Research\Research Works\My research_Time series Clustering\Data\dataset UCR\Adiac\Adiac_TEST_dismat_ED_RAW.mat';
    load(dist_mtx_file,'dismat');
    dismat=dismat+dismat';
    
    cc=hsv(k); %(length(clus{j}));
    for j=6:6%k
     mem=find(p==j);
     for i=1:5
        ind=mem(i,1);
        color=cc(j,:);
        t_traj{ind}=[1:1:length(nor_traj_raw{1})];
        plot(t_traj{ind},nor_traj_raw{ind},'color',color,'LineWidth',2)
        hold on
     end
    end
    
    %Plot_time_series(2,2,p,p,[],nor_traj_raw,[],k,5,0);
    %Plot_time_series_luminate(0,0,p,[1:1:length(p)]',[],nor_traj_raw,[],1,2,0.5,1);
end