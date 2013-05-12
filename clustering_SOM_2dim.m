function clustering_SOM_2dim()


dataset=[1,2,2,1; 2,3,1,1;1,1,2,1; 2,3,1,1;1,2,2,1; 1,3,1,1;1,2,2,1; 3,3,1,1;1,3,2,1; 2,3,1,3;1,2,2,3; 1,3,3,2;1,2,2,1; 2,3,1,2;1,1,2,2; 2,3,1,1;1,2,2,2; 1,3,2,2;2,2,2,1; 3,3,1,1;1,2,2,1; 2,2,1,3;2,2,2,3; 1,2,3,1];

% load('dataset.txt');
% dataset = dataset(2:3187, :)';
dataset = zscore(dataset);
mapWidth  = 3;
mapHeight = 2;
maxIteration = 1000;
startLearningRate = 0.1;
RandStream.setDefaultStream(RandStream('mcg16807','Seed',12575));
neurons = rand(mapWidth*mapHeight,4);
mapRadius = max (mapWidth, mapHeight) / 2;
timeConstant = maxIteration / log(mapRadius);
i = 1;
group = zeros(24,1);
for i = 1:maxIteration
    randval = randi(24);
    vector = dataset(randval,:);
    
    newmatrix = [neurons;vector];
    distance = pdist(newmatrix);
    s = squareform(distance);
    [mindist, bmu] = getMinimum(s, size(s,1));
    
    group(randval,1) = bmu;
    
    neighbourhoodRadius =  mapRadius * exp ( -i / timeConstant);
    radiusSquare = neighbourhoodRadius * neighbourhoodRadius;
    
    learningRate = startLearningRate * exp (-i / maxIteration);
    
    for n = 1:size(neurons,1)
        
        temp = pdist([neurons(bmu,:);neurons(n,:)]);
       % temp = pdist([getGridPosition(n,mapWidth); getGridPosition(bmu,mapWidth)]);
        
        distSquare = temp(1,1);
        %This means that the node location is within the radius
        if(distSquare < radiusSquare)
            
            m_dInfluence = exp(-(distSquare) / (2*radiusSquare));
            
            for c = 1:length(neurons(n,:))
                %update the weight vector
                neurons(n,c) = neurons(n,c) + m_dInfluence * learningRate * (vector(c) - neurons(n,c));
            end
        end
    end
    
end

%ASSIGN THE CLUSTERS
for i = 1:24
    
    vector = dataset(i,:);
    newmatrix = [neurons;vector];
    distance = pdist(newmatrix);
    s = squareform(distance);
    [mindist, bmu] = getMinimum(s, size(s,1));
    
    group(i,1) = bmu;
    
end

plot(neurons);
plot(group,'r.');
end

function [dis,num]= getMinimum(s,r)

[dis,num]=min(s(1:r-1,r))
end
