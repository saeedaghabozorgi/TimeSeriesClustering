
dis=zeros(5,5);
dis(1,2)=10;
dis(1,3)=10;
dis(1,4)=5;
dis(1,5)=3;
dis(2,3)=10;
dis(2,4)=10;
dis(2,5)=2;
dis(3,4)=10;
dis(3,5)=10;
dis(4,5)=10;
dis=dis+dis';
dis= squareform(dis);
Zh = linkage(dis,'average');
c = cluster(Zh,'maxclust',10);
figure;
[H,T] = dendrogram(Zh,'colorthreshold','default');
set(H,'LineWidth',2)

