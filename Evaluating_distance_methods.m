function Evaluating_distance_methods

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
    SAX_nor_traj=[];
    PAA_nor_traj=[];
    for i=1:length(nor_traj)
        SAX_nor_traj{i}= rep_SAX(nor_traj{i}, data_len, nseg, alphabet_size);
        PAA_nor_traj{i}= rep_PAA(nor_traj{i}, data_len, nseg, alphabet_size);
    end
    
    %----------------------------------------
    % Evaluate distance
    %     dis1=Mtx_Distance(nor_traj,nor_traj,'same','Euclid',alphabet_size,compression_ratio);
    %     dis2=Mtx_Distance(SAX_nor_traj,SAX_nor_traj,'same','Euclid',alphabet_size,compression_ratio);
    %     dis3=Mtx_Distance(SAX_nor_traj,SAX_nor_traj,'same','SAXminDis',alphabet_size,compression_ratio);
    %     dis4=Mtx_Distance(SAX_nor_traj,SAX_nor_traj,'same','SAXAPX',alphabet_size,compression_ratio);
    %     dis5=Mtx_Distance(PAA_nor_traj,PAA_nor_traj,'same','Euclid',alphabet_size,compression_ratio);
    %
    %     dis1=dis1/max(dis1(:));
    %     dis1=dis1+dis1';
    %     dis2=dis2/max(dis2(:));
    %     dis2=dis2+dis2';
    %     dis3=dis3/max(dis3(:));
    %     dis3=dis3+dis3';
    %     dis4=dis4/max(dis4(:));
    %     dis4=dis4+dis4';
    %     dis5=dis5/max(dis5(:));
    %     dis5=dis5+dis5';
    %
    %     dif1=(dis1-dis2).^2;
    %     differ(dataset_no,1)=sqrt(sum(dif1(:)));
    %
    %     dif2=(dis1-dis3).^2;
    %     differ(dataset_no,2)=sqrt(sum(dif2(:)));
    %
    %     dif3=(dis1-dis4).^2;
    %     differ(dataset_no,3)=sqrt(sum(dif3(:)));
    %
    %     dif4=(dis1-dis5).^2;
    %     differ(dataset_no,4)=sqrt(sum(dif4(:)));
    %
    %----------------------------------------
    %  Hierarchical
%     [c1]=do_Hierarchical_time(nor_traj,k,'Euclid','single',-1,-1,-1);
%     [c2]=do_Hierarchical_time(SAX_nor_traj,k,'Euclid','single',-1,alphabet_size,compression_ratio);
%     [c3]=do_Hierarchical_time(SAX_nor_traj,k,'SAXminDis','single',-1,alphabet_size,compression_ratio);
%     [c4]=do_Hierarchical_time(SAX_nor_traj,k,'SAXAPX','single',-1,alphabet_size,compression_ratio);
%     [c5]=do_Hierarchical_time(PAA_nor_traj,k,'Euclid','single',-1,alphabet_size,compression_ratio);
   % [c6]=do_Hierarchical_time(nor_traj,k,'DTW','single',-1,-1,-1);
    
    
    %----------------------------------------
    %  k-medoid
%     [c1,~]= do_kMediod_time (nor_traj,k,'Euclid',0,-1,-1,'RAW')  ;
%     [c2,itr]= do_kMediod_time(SAX_nor_traj,k,'Euclid',0,alphabet_size,compression_ratio,'SAX');
%     [c3,itr]= do_kMediod_time(SAX_nor_traj,k,'SAXminDis',0,alphabet_size,compression_ratio,'SAX');
%     [c4,itr]= do_kMediod_time(SAX_nor_traj,k,'SAXAPX',0,alphabet_size,compression_ratio,'SAX');
%     [c5,itr]= do_kMediod_time(PAA_nor_traj,k,'Euclid',0,alphabet_size,compression_ratio,'PAA');
%     [c6,itr]= do_kMediod_time(PAA_nor_traj,k,'DTW',0,alphabet_size,compression_ratio,'RAW');
    %----------------------------------------
    %  k-Means
    [c1,~]= do_kMeans_time (nor_traj,k,'Euclid',0,-1,-1,'RAW')  ;
    [c2,itr]= do_kMeans_time(SAX_nor_traj,k,'Euclid',0,alphabet_size,compression_ratio,'SAX');
    [c3,itr]= do_kMeans_time(SAX_nor_traj,k,'SAXminDis',0,alphabet_size,compression_ratio,'SAX');
    [c4,itr]= do_kMeans_time(SAX_nor_traj,k,'SAXAPX',0,alphabet_size,compression_ratio,'SAX');
    [c5,itr]= do_kMeans_time(PAA_nor_traj,k,'Euclid',0,alphabet_size,compression_ratio,'PAA');
%     [c6,itr]= do_kMeans_time(PAA_nor_traj,k,'DTW',0,alphabet_size,compression_ratio,'RAW');
    %----------------------------------------------------------------------
    %  k-Modes
%     [c1,~]= do_kModes_time (nor_traj,k,'Euclid',0,-1,-1,'RAW')  ;
%     [c2,itr]= do_kModes_time(SAX_nor_traj,k,'Euclid',0,alphabet_size,compression_ratio,'SAX');
%     [c3,itr]= do_kModes_time(SAX_nor_traj,k,'SAXminDis',0,alphabet_size,compression_ratio,'SAX');
%     [c4,itr]= do_kModes_time(SAX_nor_traj,k,'SAXAPX',0,alphabet_size,compression_ratio,'SAX');
%     [c5,itr]= do_kModes_time(PAA_nor_traj,k,'Euclid',0,alphabet_size,compression_ratio,'PAA');
%     [c6,itr]= do_kMeans_time(PAA_nor_traj,k,'DTW',0,alphabet_size,compression_ratio,'RAW');
    %----------------------------------------------------------------------
    % evaluate
    [~,~,~,~,~,~,~,~,~,quality]= do_Evaluate(p,c1,-1,-1,-1);
    differ(dataset_no,5)=quality;
    
    [~,~,~,~,~,~,~,~,~,quality]= do_Evaluate(p,c2,-1,-1,-1);
    differ(dataset_no,6)=quality;
    
    [~,~,~,~,~,~,~,~,~,quality]= do_Evaluate(p,c3,-1,-1,-1);
    differ(dataset_no,7)=quality;
    
    [~,~,~,~,~,~,~,~,~,quality]= do_Evaluate(p,c4,-1,-1,-1);
    differ(dataset_no,8)=quality;
    
    [~,~,~,~,~,~,~,~,~,quality]= do_Evaluate(p,c5,-1,-1,-1);
    differ(dataset_no,9)=quality;
    
    %         [~,~,~,~,~,~,~,~,~,quality]= do_Evaluate(p,c6,-1,-1,-1);
    %     differ(dataset_no,9)=quality;
    %
    %
    %     center=[];
    %     for i=1:k
    %         clus_center{i}=centre_mean(c,i,nor_traj);
    %     end
    % -- Compare with ground truth
    % Import the file
    %     p=train_data(:,1);
    %     if(min(p)==-1)
    %         p(p==-1)=2;
    %     end
    %     if(min(p)==0)
    %         p(p==0)=2;
    %     end
    %         for i=1:k
    %             cluster_center{i}=centre_mean(c1,i,nor_traj);
    %         end
    % -----------------------------
    %     p=[];
    %    [p]=do_Hierarchical_time(nor_traj,k,'Euclid','average',-1,-1,-1);
    %         for i=1:k
    %             class_center{i}=centre_mean(p,i,nor_traj);
    %         end
    %-------------------------------------
    
    %     [SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality]= do_Evaluate(p,c,nor_traj,class_center,clus_center);
    %     details(dataset_no,:)=[k,rows,data_len,SSEP,SSEC,RI,purity,BCubed,ConEntropy,f_measure,jacard,FM,quality];
    %
    
    %-----Visualization-----------------------------------
    % Plot_time_series(2,2,c,center,nor_traj,t_traj,clusterCount,1);
    % Plot_time_series(2,2,p,train_cluster_center,nor_traj,t_traj,clusterCo
    % unt,2);
    %     figure(1);
    %     t = sort(Z(:,3));
    %     th = t(size(Z,1)+2-k);
    %     [H,T,perm]= dendrogram(Z,0,'colorthreshold', th, 'orientation','right');
    %
    %     figure(3);
    %     t = sort(Z(:,3));
    %     th = t(size(Z,1)+2-k);
    %     [H,T,perm]= dendrogram(Z,0,'colorthreshold', th, 'orientation','right');
    
    %-- to plot with scroll bar
    %        figure(2);
    %        clf(2);
    %         hold off
    %         panel1 = uipanel('Parent',2);
    %         panel2 = uipanel('Parent',panel1);
    %         set(panel1,'Position',[0 0 0.95 1]);
    %         set(panel2,'Position',[0 -4 1 5]);
    %         set(gca,'Parent',panel2);
    %         s = uicontrol('Style','Slider','Parent',2,'Units','normalized','Position',[0.95 0 0.05 1],'Value',1,'Callback',{@slider_callback1,panel2});
    %
    %        cc=hsv(30);
    %        for  ind=1:30
    %         plot(t_traj{ind},nor_traj{perm(ind)}+5*ind,'color',cc(ind,:),'LineWidth',2)
    %         hold on;
    %        end
    %Plot_time_series(0,0,c1,cluster_center,nor_traj,t_traj,k,1,1);
    %Plot_time_series(0,0,p,class_center,nor_traj,t_traj,k,2,1);
    %----------------------------------------
end
end

function mmean=centre_mean(c,clusterNum,nor_traj)
t=find(c(:,1)==clusterNum);
n=length(nor_traj{1});
cluster_mem=zeros(1,n);
for j=1:size(t,1)
    cluster_mem(j,:)=nor_traj{t(j)};
end
mmean=mean(cluster_mem,1);
end

function slider_callback1(src,eventdata,arg1)
val = get(src,'Value');
set(arg1,'Position',[0 -val*5 1 5])
end