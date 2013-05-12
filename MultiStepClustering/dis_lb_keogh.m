function LB_Keogh=dis_lb_keogh(a,b,delta)
for i=1:length(a)
    if (i-delta)<1
        start=1;
    else
        start=i-delta;
    end
    if (i+delta)>length(a)
        finish=length(a);
    else
        finish=i+delta;
    end
    %Q
%     seq=a(1,start:finish);
%     Ua(1,i)=max(seq);
%     La(1,i)=min(seq);
    
    seq=b(1,start:finish);
    Ub(1,i)=max(seq);
    Lb(1,i)=min(seq);
end
LB_Keogh = sqrt(sum(sum([[a > Ub].* [a-Ub]; [a < Lb].* [Lb-a]].^2)));
end