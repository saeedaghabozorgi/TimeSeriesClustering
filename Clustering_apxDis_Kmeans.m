
function Clustering_apxDis_Kmeans()

for iter=1:100
    iter
    k=6;
    maxRow=100;
    X = rand(1,maxRow)*100;
    X=X';
    Centres=X(1:k);
    MaxIters=500;

        e = randperm(size(X,1));      % random initialization
        for i=1:k
            c(i,:)=X(e(i),:);
        end
    
    temp=zeros(maxRow,1);   % initialize as zero vector
    itr=1;
    while 1,
        itr=itr+1;
        d=dis(X,c);  % calculate objcets-centroid distances
        [z,g]=min(d,[],2);  % find group matrix g
        if g==temp,
            break;          % stop the iteration
        else
            temp=g;         % copy group matrix to temporary variable
        end
        for i=1:k
            f=find(g==i);
            if f            % only compute centroid if f is not empty
                c(i,:)=mean(X(find(g==i),:),1);
            end
        end
    end
    
    
       p=g;
       
         for i=1:k
            c(i,:)=X(e(i),:);
        end
    
    temp=zeros(maxRow,1);   % initialize as zero vector
    itr=1;
    while 1,
        itr=itr+1;
        d=Matrixmaxdis(X,c);  % calculate objcets-centroid distances
        [z,g]=min(d,[],2);  % find group matrix g
        if g==temp,
            break;          % stop the iteration
        else
            temp=g;         % copy group matrix to temporary variable
        end
        for i=1:k
            f=find(g==i);
            if f            % only compute centroid if f is not empty
                c(i,:)=mean(X(find(g==i),:),1);
            end
        end
    end
    
    
       c=g;
    
    % -------------------- visualization ------------------
       
%     hold off;
%     plot(X(1:maxRow),'g.','MarkerSize',24)
%     hold on;
%     for i=1:8
%         h = plot(find(p==i),X(p==i),'r')
%         set(h,'LineWidth',3)
%         % plot(X(p==i),'r','MarkerSize',12)
%         hold on;
%     end
%     for i=1:8
%         plot(find(c==i),X(c==i))
%         % plot(X(c==i),'MarkerSize',12)
%         hold on;
%     end
    
    %----------------------------------------------
    
    RS= calculate_Cluster_Rand_statistic(c,p);
    purity= Calculate_Cluster_Purity(c,p,2);
    
    detail(iter,1)=RS;
    detail(iter,2)=purity;
end


quality=mean(detail);
end



function Z=dis(X,Y)

for i=(1:length(X))
    for j=(1:length(Y))
        Z(i,j)=abs(X(i)-Y(j));

    end
end

end


function Z=Matrixmaxdis(X,Y)

for i=(1:length(X))
    for j=(1:length(Y))
        if floor(X(i))==floor(Y(j))
            a=[0;1];
        elseif X(i)<Y(j)
            a=[floor( X(i));ceil(Y(j))];
        else
            a=[ceil(X(i));floor(Y(j))];
        end
        Z(i,j)=pdist(a);

    end
end

end


function Y=pmaxdis(X)
t=1;
for i=(1:length(X)-1)
    for j=(i+1:length(X))
        if floor(X(i))==floor(X(j))
            a=[0;1];
        elseif X(i)<X(j)
            a=[floor( X(i));ceil(X(j))];
        else
            a=[ceil(X(i));floor(X(j))];
        end
        Y(t)=pdist(a);
         t=t+1;
    end
end

end

function Y=pmindis(X)
t=1;
for i=(1:length(X)-1)
    for j=(i+1:length(X))
        if floor(X(i))==floor(X(j))
            a=[0;0];
        elseif X(i)<X(j)
            a=[ceil( X(i));floor(X(j))];
        else
            a=[floor(X(i));ceil(X(j))];
        end
        Y(t)=pdist(a);
        t=t+1;
    end
end
end