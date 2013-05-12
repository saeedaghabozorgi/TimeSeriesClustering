function [BCubed,f_precision,f_recall]=Calculate_Cluster_BCubed_precision(c,p)
%http://fht.byu.edu/prev_workshops/workshop11/papers/1-1-schone%20Paper.pdf
correctness=zeros(length(c),length(c));
for i=1:length(c)
    for j=1:length(c)
        if(c(i)==c(j) && p(i)==p(j))
            correctness(i,j)=1;
        else
           correctness(i,j)=0; 
        end
    end
end

for i=1:length(c)
    clusterNum=c(i);
    mems=find(c==clusterNum);
    precision(i,1)=sum(correctness(i,:))/length(mems);
end
f_precision=sum(precision)/length(c);

for i=1:length(p)
    clusterNum=p(i);
    mems=find(p==clusterNum);
    recall(i,1)=sum(correctness(i,:))/length(mems);
end
f_recall=sum(recall)/length(c);
f=(1/f_precision)+(1/f_recall);
BCubed=2/f;

end