function sim=Calculate_Cluster_Similarity(c,p)
    cluster_count= length(unique(p));
    pu=[1:cluster_count]'; % cluster numbers in (p or c)
    for i=1:cluster_count
        memc=[];
        memc=find(c==i);
        pu(i,2)= length(memc);
        max_com=0;
        for j=1:cluster_count
            memp = find(p==j);
            % number of items of class j assigned to our cluster i.
            com = ismember(memp,memc);
            if sum(com)>max_com
                max_com=sum(com);
                pu(i,3)=j;
                pu(i,4)= length(memp);
                pu(i,5)=max_com; % common objects
            end
        end
    end
    sim=sum(pu(:,5))/length(c);

end

