function [C]=do_CAST_time(nor_traj_raw,dist,theroshold,varargin)
options = struct('dis',[],'alphabet_size',0,'compression_ratio',0,'rep','RAW','dis_method','Euclid');
optionNames = fieldnames(options);
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
    error('EXAMPLE needs propertyName/propertyValue pairs')
end

for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
    inpName = lower(pair{1}); %# make case insensitive
    if any(strmatch(inpName,optionNames))
        options.(inpName) = pair{2};
        %disp(['  ', num2str(inpName),' : ',num2str(options.(inpName))]);
        %    else
        %       error('%s is not a recognized parameter name',inpName)
    end
end


% representation
nor_traj=represent_TS(nor_traj_raw,options.rep,varargin{:});

C=[];
obj=nor_traj;
U=[1:length(obj)]';
C_open=[];
if ~isempty(dist)
    dis=Mtx_Distance(nor_traj,nor_traj,'same','Norm',varargin{:});
else
    dismatrix=dist
    Nor = dismatrix - min( dismatrix(:) );
    if max( Nor(:) ) ~= 0
        dismatrix = Nor / max( Nor(:) );
    else
        dismatrix=Nor;
    end
    dis=dismatrix;
end
sim=1-dis;
sim(1:length(sim)+1:length(sim)*length(sim))=0;
Cluster_num=1;
C=zeros(length(obj),1);
% only for print
Affin=sum(sim,2)./(size(sim,1)-1);
fix_val =  mean(Affin);
%disp(['       --> FIX_VAL:',num2str(fix_val),' | ',' obj:',num2str(length(obj))]);
while (~isempty(U))
    if theroshold==-6
        fix_t=theroshold6;
    elseif theroshold==-5
        fix_t=calculateT5([1:1:length(sim)]',sim);
    elseif theroshold==-4
        fix_t=calculateT4([1:1:length(sim)]',sim);
    elseif theroshold==-3
        fix_t=calculateT3(U,sim);
    elseif theroshold==-2
        fix_t=calculateT2(sim);
    elseif theroshold==-1
        fix_t=calculateT1(U,sim); % ECAST
    else
        fix_t=theroshold;
    end
    C_open=[];
    old_c=[];
    [~,inx]=MaxMat(sim,U);
    u=U(inx);
    C_open=[C_open;u];
    U(inx,:)=[];
    % Update affinity of all nodes
    a_U=sum(sim(U,C_open),2);
    a_C_open=sum(sim(C_open,C_open),2);
    iteration=0;
    while (~isequal(old_c,C_open) && iteration<10 && ~isempty(U))
        iteration=iteration+1;
        old_c=C_open;
        %addition step
        while max(a_U)>=fix_t ;
            [~,inx]= max(a_U,[],1);
            u=U(inx);
            C_open=[C_open;u];
            U(inx,:)=[];
            % Update affinity of all nodes
            a_U=sum(sim(U,C_open),2)/length(C_open);
            a_C_open=sum(sim(C_open,C_open),2)/(length(C_open)-1);
        end
        % Removal Step
        while ( min(a_C_open,[],1)<(fix_t) && min(a_C_open,[],1) >0);
            [~,inx]=min(a_C_open,[],1);
            u=C_open(inx);
            C_open(inx)=[];
            U=[U;u];
            % Update affinity of all nodes
            a_U=sum(sim(U,C_open),2)/length(C_open);
            a_C_open=sum(sim(C_open,C_open),2)/(length(C_open)-1);
        end
    end
    C(C_open)=  Cluster_num;
    Cluster_num=Cluster_num+1;
end
end


function update_affinity(U,C_open,sim)
for i=1:length(U) % Update a(U)
    a(U(i),1)=sum(sim(U(i),C_open),2);
end
for i=1:length(C_open) % Update a(C_open)
    if  u~=C_open(i)
        a(C_open(i),1)=a(C_open(i),1)+ sim(C_open(i),u);
    end
end
end


function T=calculateT1(U,sim) % ECAST
%fix_val=0.5;
Affin=sum(sim,2)./(size(sim,1)-1);
fix_val =  mean(Affin);
sim2=sim(U,U);
sim2(1:length(sim2)+1:length(sim2)*length(sim2))=0;
sim2= squareform(sim2);
sim2=sim2-fix_val;
sim2=sim2(sim2>0);
T=mean(mean(sim2))+fix_val;
if isnan(T)
    T=1;
end
end

function T=calculateT2(sim) % using the sigma
Affin=sum(sim,2)./(size(sim,1)-1);
mu = mean(Affin);
sigma = std(Affin);
outliers = (Affin - mu) < -1*sigma;
nol=sum(outliers);
if nol==0
    T=0;
else
    T=max(Affin(outliers))+0.0001;
end
end

function T=calculateT3(U,sim)
sim=sim(U,U);
Affin=sum(sim,2)./(size(sim,1)-1);
mu = mean(Affin);
sigma = std(Affin);
outliers = Affin - mu < -1*sigma;
nol=sum(outliers);
if nol==0
    T=0;
else
    T=max(Affin(outliers));
end
end

function T=calculateT4(U,main_sim) % recursive (one use or multi use)
sim=main_sim(U,U);
Affin=sum(sim,2)./(size(sim,1)-1);
mu = mean(Affin);
sigma = std(Affin);
outliers = Affin - mu < -2*sigma;
nol=sum(outliers);
if nol==0
    T=min(Affin);
else
    notOutliers=1-outliers;
    U=U(notOutliers==1);
    T=calculateT4(U,main_sim);
end
end

function T=calculateT5(U,main_sim) % single run

end



function [m,i]=MaxMat(d,U)
d=d(U,U);
d=tril(d,-1);
m=max(d);
[m,i]=max(m,[],2);
end