function [C]=do_ECAST_2dim(D,t)
D=rand(100,2);
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
% t=0.35
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


C=[];
ECAST=0;
if t==-1
    ECAST=1;
end

U=(1:length(D))';
% distance calculation
dis=dis_euclidean_matrix(D,D);

%calculate similarity
sim=1./(1+dis);
%sim=1-NorDis;

Cluster_num=1;
C=zeros(length(D),1);
while (~isempty(U))
    if ECAST==1
        t=calculateT(U,sim);
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
    for i=1:length(U) % Update a(U)
        a(U(i),1)=a(U(i),1)+ sim(U(i),u);
    end
    for i=1:length(C_open) % Update a(C_open)
        if  u~=C_open(i)
            a(C_open(i),1)=a(C_open(i),1)+ sim(C_open(i),u);
        end
    end
    iteration=0;
    while (~isequal(old_c,C_open) && iteration<200)
        iteration=iteration+1;
        old_c=C_open;
        line(D(C_open,1),D(C_open,2),'color',[.8 .1 .9],'marker','o','linestyle','none','markersize',10)
        %addition step
        while max(a(U),[],1)>=t*length(C_open);
            [~,inx]= max(a(U),[],1);
            u=U(inx);
            C_open=[C_open;u];
            hold on;
            % scatter(D(C_open,1),D(C_open,2),[],'r');
            line(D(C_open,1),D(C_open,2),'color',[.8 .1 .9],'marker','o','linestyle','none','markersize',10)
            %     pause;
            U(inx,:)=[];
            % Update affinity of all nodes
            for i=1:length(U) % Update a(U)
                a(U(i),1)=a(U(i),1)+ sim(U(i),u);
            end
            for i=1:length(C_open) % Update a(C_open)
                if  u~=C_open(i)
                    a(C_open(i),1)=a(C_open(i),1)+ sim(C_open(i),u);
                end
            end
        end
        % Removal Step
        while min(a(C_open),[],1)<t*(size(C_open,1)-1);
            [~,inx]=min(a(C_open),[],1);
            u=C_open(inx);
            line(D(u,1),D(u,2),'color',[.8 .1 .9],'marker','+','linestyle','none','markersize',10)
            C_open(inx,:)=[];
            U=[U;u];
            hold on;
            scatter(D(u,1),D(u,2),[],'b');
            %Update affinity of all nodes
            for i=1:length(U)
                if  u~=U(i)
                    a(U(i),1)=a(U(i),1)- sim(U(i),u);
                end
            end
            for i=1:length(C_open)
                a(C_open(i),1)=a(C_open(i),1)- sim(C_open(i),u);
            end
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

function [m,i]=MaxMat(d,U)
d=tril(d,-1);
d=d(:,U);
m=max(d);
[m,i]=max(m,[],2);
end

function T=calculateT(U,sim)
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