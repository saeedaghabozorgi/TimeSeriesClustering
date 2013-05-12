function clustering_kMeans_time_2Lev

k=10;
testcase=0;
all=[0,0,0,0];



%tic
data_len      = 360;
nseg          = 36;
alphabet_size = 4;
    %compression_ratio: original_data_len / symbolic_len
compression_ratio= data_len/nseg;
[nor_traj,t_traj]=ImportBankData(1,1000,'sequence',1,360);
main_nor_traj=nor_traj;



for i=1:length(nor_traj)
    SAX_nor_traj{i}= timeseries2symbol(nor_traj{i}, data_len, nseg, alphabet_size);
    %PAA_nor_traj{i}= timeseries2PAA(nor_traj{i}, data_len, nseg, alphabet_size);
end
%main_PAA_nor_traj=PAA_nor_traj;
main_SAX_nor_traj=SAX_nor_traj;
testcase=100;
while testcase<600
    testcase=testcase+1;
    nor_traj=main_nor_traj(1:testcase);
    %----------- Level 1 --------------------------------
    SAX_nor_traj=main_SAX_nor_traj(1:testcase);
   % PAA_nor_traj=main_PAA_nor_traj(1:testcase);
    
    

    
    
    
    maxRow=length(nor_traj);
    for i=1:k
        center{i}=nor_traj{i};
        center_SAX{i}=SAX_nor_traj{i};
     %   center_PAA{i}=PAA_nor_traj{i};
    end
    
    temp=zeros(maxRow,1);   % initialize as zero vector
    apitr=0;
    
    while 1,
        apitr=apitr+1;
        %     tic
        
        %dis=Mtx_SAX_APX_Distance(SAX_nor_traj,SAX_nor_traj,'same');
        %dis=Mtx_Euclid_Distance(SAX_nor_traj,center_SAX,'cell_notsame');
        %dis=Mtx_Euclid_Distance(PAA_nor_traj,center_PAA,'cell_notsame');
        %dis=Mtx_Euclid_Distance(nor_traj,center,'cell_notsame');
        %dis=Mtx_SAX_LCS_Distance(SAX_nor_traj,center_SAX,'');
        dis=Mtx_SAX_Min_Distance(SAX_nor_traj,center_SAX,'',compression_ratio);
        %t(apitr)=toc
        
        [z,c]=min(dis,[],2);  % find group matrix g
        if (c==temp | apitr==100),
            break;          % stop the iteration
        else
            temp=c;         % copy group matrix to temporary variable
        end
        for i=1:k
            center{i}=centre_mean(c,i,nor_traj,center{i});
            center_SAX{i}= timeseries2symbol(center{i}, data_len, nseg, alphabet_size);
        %    center_PAA{i}= timeseries2PAA(center{i}, data_len, nseg, alphabet_size);
        end
    end
    %----------- Level 2 --------------------------------
    itr=0;
    while 1,
        itr=itr+1;
        %dis=Mtx_SAX_Min_Distance(SAX_nor_traj,center_SAX,'');
        dis=Mtx_Euclid_Distance(nor_traj,center,'cell_not_same');
        
        [z,c]=min(dis,[],2);  % find group matrix g
        xx=c==temp;
        yy=xx==0;
        te=c(yy);
        zz=sum(yy);
        if (c==temp | itr==100),
            break;          % stop the iteration
        else
            temp=c;         % copy group matrix to temporary variable
        end
        for i=1:k
            center{i}=centre_mean(c,i,nor_traj,center{i});
        end
    end
    
    
    %---------------------------------------
    % to test if the clustering is done from scratch
    for i=1:k
        center{i}=nor_traj{i};
    end
    mitr=0;
    while 1,
        % mitr is the real iterations for  clustering from scratch
        mitr=mitr+1;
        dis=Mtx_Euclid_Distance(nor_traj,center,'cell_not_same');
        [z,c]=min(dis,[],2);  % find group matrix g
        if (c==temp | mitr==100),
            break;          % stop the iteration
        else
            temp=c;         % copy group matrix to temporary variable
        end
        for i=1:k
            center{i}=centre_mean(c,i,nor_traj,center{i});
        end
    end

    %-----------------------------------------
    
    
    
    all=[all;testcase,mitr,apitr,itr]
end
%t(1)=toc;

clusterCount=max(c)

for j=1:clusterCount
    cl_count_h(i,1)=size(find(c==i),1);
end

disp(cl_count_h);
Plot_time_series(4,3,c,center,nor_traj,t_traj,clusterCount);
end

function mmean=centre_mean(c,clusterNum,nor_traj,old_center)
t=find(c(:,1)==clusterNum);
if ~isempty (t)
    n=length(nor_traj{1});
    cluster_mem=zeros(1,n);
    for j=1:size(t,1)
        cluster_mem(j,:)=nor_traj{t(j)};
    end
    mmean=mean(cluster_mem,1);
else
    mmean=old_center;
end
end

