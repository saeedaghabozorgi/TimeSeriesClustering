function clustering_CAST_2dim_inc()
%D=rand(520,2);
%save('data\CAST_randData_2+100_inc', 'D');
load('data\CAST_randData_2+100_inc.mat', 'D');
Data=D(1:10,:);
figure(1);
hold off;
scatter(Data(:,1),Data(:,2),[],'b','c');
cc=hsv(250);
C=[];

% load('data\CAST_clusters_2+100_inc.mat', 'C');
% if ~isempty(C)
%     vv=C;
%     hold on;
%     for i=1:max(vv);
%         plot(D((vv==i),1),D((vv==i),2),'g+--');
%     end
%     for i=1:length(Data);
%         text(Data(i,1)+.01,Data(i,2), strcat('   ',num2str(i)) ,'FontSize',8,'color','black')
%     end
% end

U=(1:length(Data));
U=U';
C_open=[];
t=.8;
xx=0;
dis=dis_euclidean_matrix(Data,Data);
mindis=dis-rand()*0.1;
mm=max(max(dis),[],2);
mm=1.26;
sim=1-dis/mm;
Cluster_num=1;
C=zeros(length(Data),1);
while (~isempty(U))
    a=[];
    for i=1:length(U)
        a(U(i),1)= 0;
    end
    C_open=[];
    old_c=[];
    [m,inx]=MaxMat(sim,U);
    u=U(inx);
    C_open(end+1,1)=u;
    hold on;
    scatter(Data(C_open,1),Data(C_open,2),[],'r');
    
    %pause;
    U(inx,:)=[];
    for i=1:length(U)
        a(U(i),1)=a(U(i),1)+ sim(U(i),u);
    end
    for i=1:length(C_open)
        a(C_open(i),1)=a(C_open(i),1)+ sim(u,u);
    end
    iteration=0;
    while (~isequal(old_c,C_open) && iteration<200)
        iteration=iteration+1;
        old_c=C_open;
        %addition step
        xx=t*length(C_open);
        a_U=a(U);
        while max(a_U)>=xx
            [~,inx]= max(a_U);
            u=U(inx);
            C_open(end+1,1)=u;
            hold on;
            scatter(Data(C_open,1),Data(C_open,2),[],'r');
            %     pause;
            U(inx,:)=[];
            % Update affinity of all nodes
            for i=1:length(U)
                a(U(i),1)=a(U(i),1)+ sim(U(i),u);
            end
            for i=1:length(C_open)
                a(C_open(i),1)=a(C_open(i),1)+ sim(C_open(i),u);
            end
            xx=t*length(C_open);
            a_U=a(U);
        end
        % Removal Step
        xx=t*length(C_open);
        a_C_open=a(C_open);
        while min(a_C_open)<xx
            
            [~,inx]=min(a_C_open);
            u=C_open(inx);
            C_open(inx)=[];
            U(end+1,1)=u;
            hold on;
            scatter(Data(u,1),Data(u,2),[],'b');
            %Update affinity of all nodes
            for i=1:length(U)
                a(U(i),1)=a(U(i),1)- sim(U(i),u);
            end
            for i=1:length(C_open)
                a(C_open(i),1)=a(C_open(i),1)- sim(C_open(i),u);
            end
            xx=t*length(C_open);
            a_C_open=a(C_open);
        end
    end
    
    b(C_open)=a(C_open);
    for i=1:length(C_open)
        C(C_open(i))=  Cluster_num;
        
    end
    hold on;
    scatter(Data(C_open,1),Data(C_open,2),[],cc(Cluster_num*10,:),'filled');
    %   pause;
    Cluster_num=Cluster_num+1;
end
clusterCount=max(C);
center=zeros(clusterCount,2);
for j=1:clusterCount
    plot(Data((C==j),1),Data((C==j),2),'g+--');
    hold on;
    center(j,:)=mean(Data((C==j),:),1);
    hold on;
    scatter(center(j,1),center(j,2),[],'rs','filled');
    hold on;
    text(center(j,1)+.01,center(j,2), strcat('   ',num2str(j)) ,'FontSize',8,'color','red')
end

x=[];
U_list=[];
ddd=1;
for ii=1:300
    C_open=[];
    
    if isempty(U_list)
        x=D(10+ddd,:);
        ddd=ddd+1;
        new_node=length(Data)+1
    else
        new_node=U_list(1);
        x=Data(U_list(1),:);
        U_list(1)=[];
    end
    Data(new_node,:)=x;
    C(new_node)=0;
    b(new_node)=0;
    
    hold on;
    scatter(x(1,1),x(1,2),[],'k','filled');
    hold on;
    text(x(1,1)+.01,x(1,2), strcat('  ',num2str(new_node)) ,'FontSize',8);
    
    %---------------------------------
    % find C_open_inx by sum of similarity to a center
    %     a=zeros(1,clusterCount);
    %     dist=Mtx_Euclid_Distance(x,Data,'');
    %     sim=1-dist/mm;
    %     for j=1:clusterCount
    %         C_open= find(C==j);
    %         for i=1:length(C_open)
    %          a(j)=a(j)+sim(1,C_open(i));
    %         end
    %         a(j)= a(j)/length(C_open);
    %     end
    %     [~,C_open_inx]= max(a);
    %---------------------------------
    % find C_open_inx by nearest point

    dist=Mtx_Euclid_Distance(x,Data,'');
        unassigned_data= find(C ==0)
    dist(1,unassigned_data)=1.26;
    [~,nearest_point]= min(dist)
    C_open_inx=C(nearest_point);
    if C_open_inx==0
        rr=1;
    end
    %---------------------------------
    % find C_open_inx by nearest center
    %  dist=Mtx_Euclid_Distance(x,center,'');
    %   [~,C_open_inx]=min(dist);
    
    
    C_open= find(C==C_open_inx);
    dist=Mtx_Euclid_Distance(x,Data(C_open,:),'');
    sim=1-dist/mm;
    %add
    b(C_open)=b(C_open)+sim;
    b(new_node)=1+sum(sim);
    C(new_node)=C_open_inx;
    
    C_open= find(C==C_open_inx);
    
    a_new_clus=b(C_open);
    xx=t*length(a_new_clus);
    
    
    
    if min(a_new_clus)<xx;
        [~,inx]=min(a_new_clus);
        if C_open(inx)==new_node
            C(new_node)=Cluster_num;
            Cluster_num=Cluster_num+1;
        else
            while min(a_new_clus)<xx
                [~,inx]=min(a_new_clus);
                out_node=C_open(inx);
                x=Data(out_node,:)
                U_list=[U_list;out_node];
                C_open=find(C==C_open_inx);
                mems=Data(C_open,:);
                dist=Mtx_Euclid_Distance(x,mems,'');
                sim=1-dist/mm;
                b(C_open)=b(C_open)-sim;
                C(out_node)=0;
                C_open= find(C==C_open_inx);
                a_new_clus=b(C_open);
                xx=t*length(a_new_clus);
            end
        end
    end
    
    
    clusterCount=max(C);
    center=zeros(clusterCount,2);
    for j=1:clusterCount
        center(j,:)=mean(Data((C==j),:),1);
    end
    
    %     hold off;
    %     scatter(Data(:,1),Data(:,2),[],'b','c');
    %     if ~isempty(C)
    %         hold on;
    %         for j=1:max(C);
    %             plot(Data((C==j),1),Data((C==j),2),'g+--');
    %             hold on;
    %             scatter(Data((C==j),1),Data((C==j),2),[],cc(j*10,:),'filled');
    %             hold on;
    %             scatter(center(j,1),center(j,2),[],'rs','filled');
    %             hold on;
    %             text(center(j,1)+.01,center(j,2), strcat('  ',num2str(j)) ,'FontSize',8,'color','red');
    %         end
    %         for i=1:length(Data);
    %             text(Data(i,1)+.01,Data(i,2), strcat('  ',num2str(i)) ,'FontSize',8);
    %         end
    %     end
end

%--------------------------------------------
hold off;
scatter(Data(:,1),Data(:,2),[],'b','c');
if ~isempty(C)
    hold on;
    for j=1:max(C);
        plot(Data((C==j),1),Data((C==j),2),'g+--');
        hold on;
        scatter(Data((C==j),1),Data((C==j),2),[],cc(j*10,:),'filled');
        hold on;
        scatter(center(j,1),center(j,2),[],'rs','filled');
        hold on;
        text(center(j,1)+.01,center(j,2), strcat('  ',num2str(j)) ,'FontSize',8,'color','red');
    end
    for i=1:length(Data);
        text(Data(i,1)+.01,Data(i,2), strcat('  ',num2str(i)) ,'FontSize',8);
    end
end
%--------------------------------------------
% save('data\CAST_clusters_2+100_inc.mat', 'C');
% save('data\CAST_randData_2+100_inc.mat', 'D');
end

function [m,i]=MaxMat(d,U)
d=tril(d,-1);
d=d(:,U);
m=max(d);
[m,i]=max(m,[],2);
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