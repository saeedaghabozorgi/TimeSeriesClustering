function [c]=do_Fuzzy_time(nor_traj,k,fuzziness)
eps=.001;


maxRow=length(nor_traj);
for i=1:k
    center{i}=nor_traj{i};
end

clusterCount=k;
FClusDisMatrix=Mtx_Euclid_Distance(nor_traj,center,'cell_not_same');
OffMembership=zeros(length(nor_traj),clusterCount);

for i=1:length(nor_traj)
    denominator=0;
    for j=1:clusterCount
        if FClusDisMatrix(i,j)==0 ;
            FClusDisMatrix(i,j)=eps;
        end
        denominator=denominator+(1/FClusDisMatrix(i,j))^fuzziness;
    end
    for j=1:clusterCount
        mem=((1/FClusDisMatrix(i,j))^fuzziness)/denominator;
        f0(i,j)=mem;
    end
end



N=length(nor_traj);
f = zeros(N,clusterCount);                % partition matrix
iter = 0;                       % iteration counter

while  max(max(f0-f)) > .01
    iter = iter + 1;
    f = f0;
    
    % recalculate the centers
    fm = f.^fuzziness;
    sumf = sum(fm);
    %  v = (fm'*data.X)./(sumf'*ones(1,2));
    for j=1:clusterCount
        numerator=zeros(1,length(nor_traj{1}));
        for i=1:length(nor_traj)
            numerator=numerator+(fm(i,j)*nor_traj{i});
        end
        %numerator=(fm(:,j)'*data.X);
        denominator=sumf(1,j);
        center{j}=numerator/denominator;
    end
    
    %update dist
    FClusDisMatrix=Mtx_Euclid_Distance(nor_traj,center,'cell_not_same');
    
    
    % Update f0
    d = (FClusDisMatrix+1e-10).^(-1/(fuzziness-1));
    f0 = (d ./ (sum(d,2)*ones(1,clusterCount)));
    
    % re assigning to clusteres
    [~,c]=max(f,[],2);
    
    

    
    % disp(cl_count_h);
    % Plot_time_series(5,4,c,center,nor_traj,t_traj,clusterCount);
    %   pause;
    
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


