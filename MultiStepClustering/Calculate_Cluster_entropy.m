function entropy=Calculate_Cluster_entropy(c,p)

%Demo

% c=[1 ;1; 2; 2; 2; 3; 3; 4];
% p=[1; 2 ;2; 2; 3; 4 ;3 ;2];

n=size(c,1); % size of cluster

%calculate the entropy of each cluster (Ej)
cluster_count=length(unique(c));
for j=1:cluster_count
    H=0;
    mem_c= find(c==j); % Mems of Cj
    class_count=length(unique(p));
    % for each class
    for i=1:class_count
        % calculate probability Pr(Pi|Cj), that an instance in Cj belongs to
        % class	Pi , where Pr(Pi|Cj) = |Cj?Pi|/|Cj|.
        mem_p=find(p==i);
        Pr=sum(ismember(mem_c,mem_p));
        Pr=Pr/size(mem_c,1);
        if Pr~=0  % TO DO: if X and Y are independent, then knowing X does not give any information about Y and vice versa, so their mutual information is zero.
        H=H+Pr*log(Pr);
        end
    end
    E(j)= -(1/log(class_count))*H;  % Ej=normalized entropy Cj
end
% The overall entropy (E ? [0, 1]) is defined as the sum of the individual cluster entropies weighted by the size of each cluster:

pE=0;
for j=1:cluster_count
    mem_c= find(c==j); % Mems of Cj
    pE=pE+ size(mem_c,1)*E(j);
end
entropy= pE/size(c,1);
end

