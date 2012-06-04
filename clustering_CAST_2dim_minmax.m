function clustering_CAST_2dim_minmax()
load('data\data_minmax.mat','D');
%D=[.1 .2; .5 .9; .2 .3; .8 .8; .2 .45];
%D=rand(50,2);
load('data\clusters_minmax.mat', 'C');

vv=C;
figure(1);
hold off;
scatter(D(:,1),D(:,2),[],'b','c');

cc=hsv(200);

  hold on;
  for i=1:max(vv);
    plot(D((vv==i),1),D((vv==i),2),'g+--');
  end
for i=1:length(D);
  text(D(i,1)+.02,D(i,2), strcat('   ',num2str(i)) ,'FontSize',10)
end
%D=[1;0;2;3; 7;9;8 ;10;11;13];
U=(1:length(D));
U=U';
C_open=[];
t=.75;
xx=0;

dis=dis_euclidean_matrix(D,D);
mindis=squareform(dis);
mindis=mindis -rand(1,length(mindis))*.300;
mindis=squareform(mindis);


maxdis=squareform(dis);
maxdis=maxdis +rand(1,length(maxdis))*.300;
maxdis=squareform(maxdis);

mm=1.26;
nmm=max(max(mindis),[],2);
xmm=max(max(maxdis),[],2);

 sim=1-dis/mm;
 maxsim=1-mindis/nmm;
 minsim=1-maxdis/xmm;

% sim=(maxsim+minsim)/2;
 sim=maxsim;

Cluster_num=1;
C=zeros(length(D),1);
while (~isempty(U))
      aa=[];
    for i=1:length(U)
        aa(U(i),1)= 0;
    end
    C_open=[];
    old_c=[];
    [m,inx]=MaxMat(sim,U);
   % [m,inx]=MaxMat(minsim,U);
    u=U(inx);
    C_open(end+1,1)=u;
    hold on;
    scatter(D(C_open,1),D(C_open,2),[],'r');
 
    %pause;
    U(inx,:)=[];
    for i=1:length(U)
       aa(U(i),1)=aa(U(i),1)+ sim(U(i),u);
         %aa(U(i),1)=aa(U(i),1)+ minsim(U(i),u);
    end
    for i=1:length(C_open)
        %aa(C_open(i),1)=aa(C_open(i),1)+ minsim(C_open(i),u);
        aa(C_open(i),1)=aa(C_open(i),1)+ sim(C_open(i),u);
    end
    iteration=0;
    while (~isequal(old_c,C_open) && iteration<200)
        iteration=iteration+1;
        old_c=C_open;
        %addition step
        xx=t*length(C_open);
        a_U=aa(U);
        while max(a_U)>=xx
            [~,inx]= max(a_U);
            u=U(inx);
            C_open(end+1,1)=u;
            hold on;
            scatter(D(C_open,1),D(C_open,2),[],'r');
        %     pause;
            U(inx,:)=[];
            % Update affinity of all nodes
            for i=1:length(U)
               % aa(U(i),1)=aa(U(i),1)+ minsim(U(i),u);
                aa(U(i),1)=aa(U(i),1)+ sim(U(i),u);
            end
            for i=1:length(C_open)
               % aa(C_open(i),1)=aa(C_open(i),1)+ minsim(C_open(i),u);
                aa(C_open(i),1)=aa(C_open(i),1)+ sim(C_open(i),u);
            end
            xx=t*length(C_open);
            a_U=aa(U);
        end
        % Removal Step
        xx=t*length(C_open);
        a_C_open=aa(C_open);
        while min(a_C_open)<xx
            [~,inx]=min(a_C_open);
            u=C_open(inx);
            C_open(inx)=[];
            hold on;
            scatter(D(u,1),D(u,2),[],'b');
            U(end+1,1)=u;
            %Update affinity of all nodes
            for i=1:length(U)
               % aa(U(i),1)=aa(U(i),1)- maxsim(U(i),u);
                aa(U(i),1)=aa(U(i),1)- sim(U(i),u);
            end
            for i=1:length(C_open)
               % aa(C_open(i),1)=aa(C_open(i),1)- maxsim(C_open(i),u);
                aa(C_open(i),1)=aa(C_open(i),1)- sim(C_open(i),u);
            end
            xx=t*length(C_open);
            a_C_open=aa(C_open);
        end
    end
    
    
    for i=1:length(C_open)
        C(C_open(i))=  Cluster_num;
    end
     hold on;
     scatter(D(C_open,1),D(C_open,2),[],cc(Cluster_num*10,:),'filled');
    % plot(D(C_open,1),D(C_open,2),'g+--')

  %   pause;
    Cluster_num=Cluster_num+1;
end
%save('data\clusters_minmax.mat', 'C');
%save('data\data_minmax.mat', 'D');
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