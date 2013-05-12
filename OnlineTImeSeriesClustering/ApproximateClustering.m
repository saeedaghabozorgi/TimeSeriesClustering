function test
% to create the clusters of one time series
TS = load('Bird_MFCC.txt');
tL = int64(50000);
mL = int64(50);
wL = int64(500);

T = double(tL);
M = double(mL);
W = double(wL);


scrsz = get(0,'ScreenSize');


h = subplot(5,5,1:4);
% h1 = subplot(4,2,3);
% h2 = subplot(4,2,4);
% h3 = subplot(4,2,5);
% h4 = subplot(4,2,6);
% h5 = subplot(4,2,7);
% h6 = subplot(4,2,8);

Mark = zeros(1,length(TS));
i = 501;
s=1;
dis=inf(200,200);
clsCount=20

Cr=[1:1:clsCount];
buk=[];
pos=1;
edd=[0];
while pos<length(TS)-M
    
    if isempty(buk)
        buk(1).elm(1,:)=[pos M 0];
        buk(1).cen= DNorm_Unif(TS(1:M),8);
    else
        bsf_ed=10;
        for bi=1:length(buk)
            ed = ED(buk(bi).cen,DNorm_Unif(TS(pos:pos+M-1),8));
            edd(1,end+1)=ed;
            if (bsf_ed > ed)
                bsf_ed  = ed;
                bsf_pos = pos;
                bsf_bi = bi;
            end
        end
        if bsf_ed<10
            buk(bsf_bi).elm(end+1,:) = [bsf_pos M bsf_ed];
        else
            blen = length(buk);
            buk(blen+1).elm(1,:)=[pos M 0];
            buk(blen+1).cen= DNorm_Unif(TS(pos:pos+M-1),8);
        end
    end


pos=pos+5;

end
for bi=1:length(buk)
    a(bi,1)=size(buk(bi).elm,1);
    a(bi,2)=mean(buk(bi).elm(:,3));
end
plot(a);
end
%% Note: motif is a struct
function [cost_create]= CostIfCreate(TS,Mark,Cluster,pos)

clen = length(Cluster);
a = pos;
M = 50;
ts1 = DNorm_Unif(TS(a:a+M-1));
%%% BUG_ROUND
center = round((ts1+ts1)/2);

Mark(a:a+M-1)=1;


Xshf=0;
Cluster(clen+1).elm(1,:) = [a, M, Xshf];

Cluster(clen+1).cenTS = center;
Cluster(clen+1).cenM = M;
Cluster(clen+1).minXs = 0;
Cluster(clen+1).maxXs = 0;
Cluster(clen+1).minM = M;
Cluster(clen+1).maxM = M;

%%% New MDL: Bits Save / length
%     diff1 = MDL(ts1-center);
%     diff2 = MDL(ts2-center);
%     MDL_OLD = (MDL(ts1)+MDL(ts2))/M;
%     MDL_NEW = (MDL(center)+ min(diff1,diff2))/M;  % ????.
MDL_OLD = MDL(DNorm_Unif(TS(a:a+M-1)))/M;   %DL(A)
cost_create = MDL_OLD;
end
%%
function [Cluster Mark] = CreateNewCluster(TS,Mark,Cluster,SS)

clen = length(Cluster);
a = SS(1,2);
b = SS(1,3);
M = SS(1,4);
ts1 = DNorm_Unif(TS(a:a+M-1));
ts2 = DNorm_Unif(TS(b:b+M-1));
center = (ts1+ts2)/2;

Xshf=0;
Cluster(clen+1).elm(1,:) = [a, M, Xshf];
%Cluster(clen+1).elm(2,:) = [b, M, Xshf];
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
%%
function [Cost_Add] = CostIfAddML(TS,Mark,Cluster,Cost_Add,MM,pos)
if (length(Cluster) < 1),
    Cost_Add = [inf 0 0 0];
    return;
end

for ci=1:length(Cluster)
    [Cost_Add(ci,:)] = FindAddedClusterML(TS,Mark,Cluster,ci,MM,pos);
    cadd = Cost_Add(ci,:);
    fprintf('Add=> cid:%d, ts:%d, M=%d,  cost=%.2f\n',cadd(3),cadd(2),cadd(4),cadd(1));
end
end

%% Brute Force to find the next candidate to merge
function [cost_Add] = FindAddedClusterML(TS,Mark,Cluster,ci,MM,pos)

cenTS = Cluster(ci).cenTS;
M = Cluster(ci).cenM;

bsf_ed = inf;

ed = ED(cenTS,DNorm_Unif(TS(pos:pos+M-1)));
if (bsf_ed > ed)
    bsf_ed  = ed;
    bsf_pos = pos;
    bsf_cid = ci;
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
new_cen = ( (n-1)*old_cen + DNorm_Unif(TS(pos:pos+M-1)) )/n;
Cluster(cid).cenTS = new_cen;
end




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%% MDL Stuffs  %%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function cost = MDL_Cluster_i(TS,clusi)
cost = 0;
cenTS = clusi.cenTS;
cenM  = clusi.cenM;


cost_diff = 0;
max_mdl_diff = -inf;

for j=1:size(clusi.elm,1) %%%% Careful Here %%%
    pos  = clusi.elm(j,1);
    M    = clusi.elm(j,2);
    Xshf = clusi.elm(j,3);
    
    ts = DNorm_Unif(TS(pos:pos+M-1));
    ts = [zeros(1,Xshf), ts, zeros(1,cenM-M-Xshf)];
    
    mdl_diff = MDL(ts-cenTS);
    max_mdl_diff = max(max_mdl_diff, mdl_diff);
    cost_diff = cost_diff + mdl_diff;
end
%cost_diff = cost_diff - max_mdl_diff;
%cost=MDL(cenTS)
% cost=cost_clus+cost_diff
% DLC(C)=DL(H)+z[DL(A|H)]-max(DL(A|H)
cost = MDL(cenTS) + cost_diff - max_mdl_diff;
end

%%
function cost = MDL(ts)
cost = bitcost_ent(ts);
end

function cost = bitcost_ent(ts)
ts = ts-min(ts)+1;
if (min(ts)==max(ts))
    % ts(end) = -ts(end);
    cost=0;
else
    cost = ent(ts);
end
end

function cost = ent(ts)
min_val = min(ts);
max_val = max(ts);
f = histc(ts(:),min_val:max_val);
f = f(:)'/sum(f);
ent = -f.*log(f+0.000000000000000000001)./log(2);
cost = length(ts)* sum(ent);
end
