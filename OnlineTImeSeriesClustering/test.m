function test

d = load('LSF5_10.txt');
tL = int64(50000);
mL = int64(50);
wL = int64(200);

T = double(tL);
M = double(mL);
W = double(wL);


scrsz = get(0,'ScreenSize');
figure;


h = subplot(5,5,1:4);
% h1 = subplot(4,2,3);
% h2 = subplot(4,2,4);
% h3 = subplot(4,2,5);
% h4 = subplot(4,2,6);
% h5 = subplot(4,2,7);
% h6 = subplot(4,2,8);

Mark = zeros(1,length(d));
i = 700;
s=1;
dis=inf(200,200);
clsCount=20

Cr=[1:1:clsCount];
s=1;
Cluster=[];

for s=1:clsCount
    %  CC(s,:)=zNorm(d(s:(s+M-1)));
    CC(s,:) = [inf s s 50];
    [Cluster Mark] = CreateNewCluster(d,Mark,Cluster,CC(s,:));
end

for cid=1:clsCount
    h1 = subplot(5,5,cid+4);
    plot(h1,Cluster(cid).cenTS,'r','LineWidth',2);
end

s=1;

while i<5000
    
    plot(h,(i-W:i),d(i-W:i));
    xlim(h,[i-W,i]);
    
    hold (h, 'on');
    l1=i-W;
    plot(h,(l1:l1+M-1),d(l1:(l1+M-1)),'r','LineWidth',2);
    Sub(s,:)=zNorm(d(l1:(l1+M-1)));
    for cid=1:clsCount
        dis(s,cid)=ED(Cluster(cid).cenTS,Sub(s,:));
    end
    [dist,clusloc]=min(dis(s,:));
    SS(s,:) = [inf clusloc l1 50];
    %     [Cluster Mark] = CreateNewCluster(d,Mark,Cluster,SS);
    pos=l1;
%     
%     [Cluster Mark] = AddToCluster(d,Mark,Cluster,pos,clusloc,M);
%     s=s+1;
%     
%     h1 = subplot(5,5,clusloc+4);
%     hold off;
%     for j=1:size(Cluster(clusloc).elm,1)
%         pos=Cluster(clusloc).elm(j,1);
%         r=zNorm(d(pos:pos+M-1));
%         plot(h1,r,'b','LineWidth',2);
%         hold on;
%     end
%     plot(h1,Cluster(clusloc).cenTS,'r','LineWidth',2);
        drawnow
    i=i+5;
end



end
%%
function [Cost_Add] = CostIfAddML(TS,Mark,Cluster,Cost_Add,MM)           
    if (length(Cluster) < 1),
        Cost_Add = [inf 0 0 0];
        return;
    end
    
    for ci=1:length(Cluster) 
        [Cost_Add(ci,:)] = FindAddedClusterML(TS,Mark,Cluster,ci,MM);
        cadd = Cost_Add(ci,:);
        fprintf('Add=> cid:%d, ts:%d, M=%d,  cost=%.2f\n',cadd(3),cadd(2),cadd(4),cadd(1));         
    end
end

%% Brute Force to find the next candidate to merge
function [cost_Add] = FindAddedClusterML(TS,Mark,Cluster,ci,MM)     
    
    cenTS = Cluster(ci).cenTS;       
    M = Cluster(ci).cenM;
    
    bsf_ed = inf;
    for p=1:length(TS)-M+1
        ChkM = [0,  MM(MM<=M)-1];
        if (any(Mark(p+ChkM)==1))
            continue;
        end      
        

        ed = ED(cenTS,DNorm_Unif(TS(p:p+M-1)));
        if (bsf_ed > ed)
            bsf_ed  = ed;                
            bsf_pos = p;
            bsf_cid = ci; 
        end                
    end

    if (bsf_ed==inf)
       cost_Add = [inf 0 0 0];
    else
        %%% Old MDL: Bits Save / length
        MDL_OLD1 = MDL_Cluster_i(TS,Cluster(bsf_cid))/M;   %DLC(C)
        MDL_OLD2 = MDL(DNorm_Unif(TS(bsf_pos:bsf_pos+M-1)))/M;   %DL(A)
        MDL_OLD = MDL_OLD1 + MDL_OLD2;        % DL(Before)
        
        cen1 = Cluster(bsf_cid).cenTS;
        
        Cluster = AddToCluster(TS,Mark,Cluster,bsf_pos,bsf_cid,M);
        cen2 = Cluster(bsf_cid).cenTS;
        
        %%% New MDL: Bits Save / length
       MDL_NEW = MDL_Cluster_i(TS,Cluster(bsf_cid))/M;
       cost = MDL_NEW - MDL_OLD;              
       cost_Add = [cost bsf_pos bsf_cid M];
    end
end


%%
function [Cluster Mark] = AddToCluster(TS,Mark,Cluster,pos,cid,M)
Mark(pos:pos+M-1)=1;
Cluster(cid).elm(end+1,:) = [pos M 0];
n = size(Cluster(cid).elm,1);
old_cen = Cluster(cid).cenTS;
new_cen = ( (n-1)*old_cen + zNorm(TS(pos:pos+M-1)) )/n;
Cluster(cid).cenTS = new_cen;
end


%%
function [Cluster Mark] = CreateNewCluster(TS,Mark,Cluster,SS)

clen = length(Cluster);
a = SS(1,2);
b = SS(1,3);
M = SS(1,4);
ts1 = zNorm(TS(a:a+M-1));
ts2 = zNorm(TS(b:b+M-1));
center = (ts1+ts2)/2;

Xshf=0;
Cluster(clen+1).elm(1,:) = [a, M, Xshf];
Cluster(clen+1).elm(2,:) = [b, M, Xshf];
Cluster(clen+1).cenTS = center;
Cluster(clen+1).cenM = M;
Cluster(clen+1).minXs = 0;
Cluster(clen+1).maxXs = 0;
Cluster(clen+1).minM = M;
Cluster(clen+1).maxM = M;

%     if (any(Mark(a:a+M-1)) || any(Mark(b:b+M-1)))
%         fprintf('Warning!! CreateNewCluster:: Mark twice!!\n');
%         pause(1);
%     end
%
%     Mark(a:a+M-1)=1;
%     Mark(b:b+M-1)=1;
end


function d=ED(ts1,ts2)
d = sqrt(sum((ts1-ts2).^2));
end