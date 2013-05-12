function SSE=Calculate_Cluster_SSE(c,nor_traj,center)

% calculate SSE
clusterCount=max(c);
SSE=0;
for j=1:clusterCount
    mem=find(c==j);
    for i=1: length(mem)
        d=dis_euclidean(nor_traj{mem(i)},center{j});
        SSE=SSE+d^2;
    end
end
end