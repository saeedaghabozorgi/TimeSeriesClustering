function CSM=Calculate_Cluster_CSM(c,p)
% To compute purity , each cluster is assigned to the class which is most
% frequent in the cluster, and then the accuracy of this assignment is
% measured by counting the number of correctly assigned documents and
% dividing by length(c)

% High purity is easy to achieve when the number of clusters is large - in
% particular, purity is 1 if each document gets its own cluster. Thus,
% we cannot use purity to trade off the quality of the clustering against
% the number of clusters.



% based on http://www.cse.iitm.ac.in/~cs672/purity.pdf

% the result and approach are exactly same in both approaches


cluster_count= length(unique(p));
pu=[1:cluster_count]'; % cluster numbers in p
for i=1:cluster_count
    memc=[];
    memc=find(c==i);
    sim=0;
    max_sim=0;
    for j=1:cluster_count
        memp = find(p==j);
        % number of items of class j assigned to our cluster i.
        com = ismember(memp,memc);
        sim=2*sum(com)/(length(memc)+length(memp));
        if sim>max_sim
            max_sim=sim;
            pu(i,3)=sim; % 
            pu(i,2)=j;
        end
    end
end
CSM=sum(pu(:,3))/max(c);
end


