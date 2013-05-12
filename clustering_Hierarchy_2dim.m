function clustering_Hierarchy_2dim()
Points=rand(20,2);

figure(1);
hold off;
scatter(Points(:,1),Points(:,2),[],'b','c');
cc=hsv(100);
C=[];

k=5;

d= pdist(Points);
Z = linkage(d);
Y = inconsistent(Z);
Y = inconsistent(Z,2);
Y = inconsistent(Z,3);
%c = cluster(Z,'maxclust',k)
c = cluster(Z,.7)
figure(2);
hold off;
t = sort(Z(:,3));
th = t(size(Z,1)+2-k);
%dendrogram(Z,0,'colorthreshold', th);
dendrogram(Z,0);

figure(1);
      hold on;
   vv=c;
      for i=1:max(vv);
        plot(Points((vv==i),1),Points((vv==i),2),'g+--');
    end

end
