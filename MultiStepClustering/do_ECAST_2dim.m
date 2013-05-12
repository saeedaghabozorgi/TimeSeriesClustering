function [C]=do_ECAST_2dim(D,t)
D=rand(10,2);
t=0.7;
%save('..\data\CAST_randData_10_1.mat', 'D');
%load('..\data\CAST_randData_10_1.mat');

% why we need to calculate exact distance for all of points wheras some
% distances never used, maybe it is adequate to just calculate the exact
% distance for some and approximate

%min_aff_clus defines the minimum affilation of each cluster after clustering;
%by this parameter we can predefine the number of clusters.

%----------demo 1---- By click
%-- load
% load('..\data\CAST_randData_10_3.mat');
% clf;
% scatter(D(:,1),D(:,2),[],'b','c');
% t=0.5
%-- Save
% D=[];
% hFig = figure(1);
% axis([0 5 0 5])
% hold off
% for i=1:28
%     [x,y] = ginput(1);
%     D=[D;x,y];
%     hold on;
%     scatter(D(:,1),D(:,2),[],'b','c');
% end
%  save('..\data\CAST_randData_10_3.mat', 'D');
%-- manual
% D=[0 0;1 1; 2 1; 0 4; 2 5; 3 5; 3 4; 4 5; 4 4];
% clf;
% scatter(D(:,1),D(:,2),[],'b','c');
%----------------

cc=hsv(100);
min_aff_clus=200;
figure(1);
hold off;
scatter(D(:,1),D(:,2),[],'b','c');
cc=hsv(100);

%------
%find T based on KNN
% [neighbors distances] = kNearestNeighbors(D, D, 10);
% distances=distances(:,2:end);
% ssm=1./(1+distances);
% t=mean(mean(ssm));
theroshold=-1

C=[];


U=(1:length(D))';
% distance calculation
dis=dis_euclidean_matrix(D,D);

%calculate similarity
sim=1./(1+dis);
sim(1:length(sim)+1:length(sim)*length(sim))=0;

Cluster_num=1;
C=zeros(length(D),1);
while (~isempty(U))
    if theroshold==-5
        fix_t=calculateT4(U,sim);
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
    a=zeros(size(D,1),1);
    C_open=[];
    old_c=[];
    [~,inx]=MaxMat(sim,U);
    u=U(inx);
    C_open=[C_open;u];
    hold on;
    %scatter(D(C_open,1),D(C_open,2),[],'r');
    line(D(C_open,1),D(C_open,2),'color',[.5 .1 .9],'marker','p',  'linestyle','none','markersize',20)
    %pause;
    U(inx,:)=[];
    a_U=sum(sim(U,C_open),2);
    a_C_open=sum(sim(C_open,C_open),2);
    iteration=0;
    while (~isequal(old_c,C_open) && iteration<200 && ~isempty(U))
        iteration=iteration+1;
        old_c=C_open;
        line(D(C_open,1),D(C_open,2),'color',[.8 .1 .9],'marker','o','linestyle','none','markersize',10)
        %addition step
        while max(a_U)>=fix_t ;
            [~,inx]= max(a_U,[],1);
            u=U(inx);
            C_open=[C_open;u];
            U(inx,:)=[];
            hold on;
            % scatter(D(C_open,1),D(C_open,2),[],'r');
            line(D(C_open,1),D(C_open,2),'color',[.8 .1 .9],'marker','o','linestyle','none','markersize',10)
            % Update affinity of all nodes
            a_U=sum(sim(U,C_open),2)/length(C_open);
            a_C_open=sum(sim(C_open,C_open),2)/(length(C_open)-1);
        end
        % Removal Step
        while min(a_C_open,[],1)<(fix_t);
 [~,inx]=min(a_C_open,[],1);
            u=C_open(inx);
            C_open(inx)=[];
            U=[U;u];
            line(D(u,1),D(u,2),'color',[.8 .1 .9],'marker','+','linestyle','none','markersize',10)
            hold on;
            scatter(D(u,1),D(u,2),[],'b');
            %Update affinity of all nodes
            a_U=sum(sim(U,C_open),2)/length(C_open);
            a_C_open=sum(sim(C_open,C_open),2)/(length(C_open)-1);
        end
    end
    
    for i=1:length(C_open)
        C(C_open(i))=  Cluster_num;
    end
    %--- figure
    hold on;
    scatter(D(C_open,1),D(C_open,2),[],cc(Cluster_num*5,:),'filled');
    %     %   pause;
    %     vv=C;
    %     for i=1:max(vv);
    %         plot(D((vv==i),1),D((vv==i),2),'g+--');
    %     end
    Cluster_num=Cluster_num+1;
end
%save('data\CAST_clusters_10_1.mat', 'C');
end


function T=calculateT1(U,sim) % mycast
%fix_val=0.5;
Affin=sum(sim,2)./(size(sim,1)-1);
ma = mean(Affin);
ms=mean(squareform(sim));
fix_val = ma;
sim2=sim(U,U);
ms=mean(squareform(sim2));
sim2(1:length(sim2)+1:length(sim2)*length(sim2))=0;
sim2= squareform(sim2);
sim2=sim2-fix_val;
sim2=sim2(sim2>0);
T=mean(mean(sim2))+fix_val;
if isnan(T)
    T=0;
end
end

function [m,i]=MaxMat(d,U)
d=tril(d,-1);
d=d(:,U);
m=max(d);
[m,i]=max(m,[],2);
end

function T=calculateT(U,sim) %ECAST
ea = 0;
ecount = 0;

for i=1:length(U)-1
    for j=i+1:length(U)
        ui=U(i);
        uj=U(j);
        a(i)=sim(ui,uj);
        if a(i)>= 0.5
            ea=ea+(a(i)-0.5);
            ecount=ecount+1;
        end
    end
end
if ecount >0
    T = ( ea/ecount ) + 0.5;
else
    T=1;
end

end

function ax=affinity(x,C_open,sim)
dd=0;
for i=1:length(C_open)
    dd=dd+sim(C_open(i),x);
end
if (~isempty(C_open))
    ax=dd/length(C_open);
else
    ax=0;
end

end