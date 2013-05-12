function [C]=do_ECAST_2dim(D,t)
%D=rand(100,2);
%t=-1;
%save('..\data\CAST_randData_10_1.mat', 'D');
%load('..\data\CAST_randData_10_1.mat');

% why we need to calculate exact distance for all of points wheras some
% distances never used, maybe it is adequate to just calculate the exact
% distance for some and approximate

%min_aff_clus defines the minimum affilation of each cluster after clustering;
%by this parameter we can predefine the number of clusters.

%----------demo 1---- By click
%-- load
load('..\data\CAST_randData_10_1.mat');
clf;
scatter(D(:,1),D(:,2),[],'b','c');
t=0.92
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
% figure(1);
% hold off;
% scatter(D(:,1),D(:,2),[],'b','c');
% cc=hsv(100);

%------
%find T based on KNN
% [neighbors distances] = kNearestNeighbors(D, D, 10);
% distances=distances(:,2:end);
% ssm=1./(1+distances);
% t=mean(mean(ssm));

rrr=[];
C=[];
ECAST=0;
if t==-1
    ECAST=1;
end

U=(1:length(D))';

C_open=[];

xx=0;
% distance calculation
dis=dis_euclidean_matrix(D,D);
sim=1./(1+dis);

Cluster_num=1;
C=zeros(length(D),1);
while (~isempty(U))
    if ECAST==1
        t=calculateT(U,sim);
    end
    % TO DO  if T==Nan id U=2
    a=zeros(size(D,1),1);
    C_open=[];
    old_c=[];
    [m,inx]=MaxMat(sim,U);
    u=U(inx);
    C_open=[C_open;u];
    hold on;
    %scatter(D(C_open,1),D(C_open,2),[],'r');
    line(D(C_open,1),D(C_open,2),'color',[.5 .1 .9],'marker','p',  'linestyle','none','markersize',20)
    %pause;
    U(inx,:)=[];
    a=update_affinity (a,u,C_open,sim,5,1)
    iteration=0;
    a_U=a(U);
    [val inx]=max(a_U,[],1);
    oldVal=val;
    while abs(val-oldVal)<.06
        oldVal=val;
        rrr=[rrr;val];
        u=U(inx);
        C_open=[C_open;u];
        hold on;
        % scatter(D(C_open,1),D(C_open,2),[],'r');
        line(D(C_open,1),D(C_open,2),'color',[.8 .1 .9],'marker','o','linestyle','none','markersize',10)
        U(inx,:)=[];
        %  a=update_affinity (a,u,C_open,sim,round(size(C_open,1)/2),1)
        a=update_affinity (a,u,C_open,sim,10,1)
        a_U=a(U);
        [val inx]=max(a_U,[],1);
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

function a=update_affinity (a,u,c_open,sim,k,factor)
if(size(c_open,1)<k)
    k=size(c_open,1);
end
for i=1:size(a,1)
    simi= sim(c_open,i);
    [sortval sortpos] = sort(simi,'descend');
    neighborIds(i,:) = sortpos(1:k);
    neighborSimilarity(i,:) =sortval(1:k);
    a(i,1)= mean(sortval(1:k));
end
end

