function [labels, cost, medoids] = do_kMedoids_keogh(k,Dist)
%Dist is the output from Pdist, a vector

D = squareform(Dist);
n = size(D,1);

medoids = randsample(n,k);
[costs,labels] = min(D(medoids,:));
cost = sum(costs);

last = 0;

if length(medoids)~=k
    disp('set the correct medoids!');
    return;
end


while any(last ~= medoids)
    best_so_far_medoids = medoids;
    for i = 1:k %for each medoids point
        medoids_aux = medoids;
        for j = 1:n %for each non-medoids point
            if ismember(j, medoids)
                continue
            end
            medoids_aux(i) = j; %swap
            [costs_aux,labels_aux] = min(D(medoids_aux,:));
            cost_aux = sum(costs_aux);
            if (cost_aux < cost)
                best_so_far_medoids = medoids_aux;
                cost = cost_aux;
                labels = labels_aux;
            end
        end
    end
    last = medoids;
    medoids = best_so_far_medoids;
end

labels = labels';
