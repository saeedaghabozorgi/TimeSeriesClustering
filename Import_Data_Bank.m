function [nor_traj,t_traj]=Import_Data_Bank(FromUser,ToUser,cond,SinceTime,ToTime)
%cond: full     : from 0 to 365
%cond: time     : from TSstart to TSend
%cond: timeSeq  : from TSstart to ToTime
%cond: sequence : from SinceTime to ToTime
% Import the file
fileName='..\dataset Bank\data10.txt';

walking=1;
lng=365;
newData1 = importdata(fileName);
% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
data=newData1.(vars{1});



org_traj=[];  % orginal traj
nor_traj=[]; % normalized orginal traj
t_traj=[]




for z=FromUser:ToUser
    CurTS=data(data(:,1)==z,:);
    Obs=size(CurTS,1);   %number of observations/transacctions
    TSstart=CurTS(1,4);  %Actual Start point of trajectory
    TSend=CurTS(Obs,4);  %Actual End point of trajectory
    a=zeros(1,lng);
    % a(1:TSstart)=CurTS(1,6);
    a(1:TSstart)=0;
    
    
    
    if walking==1  % if it is based on slow walking and slops
        for i=1:Obs-1;
            tx1=CurTS(i,4);
            y1=CurTS(i,6);
            a(tx1)=y1;
            if(length(CurTS)>1)
                tx2=CurTS(i+1,4);
                y2=CurTS(i+1,6);
                a(tx2)=y2;
                dif=tx2-tx1-1;
                if dif>0
                    for j=1:1:dif
                        slope=(y2-y1)/(tx2-tx1);
                        bb=y1-slope*tx1;
                        y=slope*(tx1+j)+bb;
                        a(tx1+j)=y;
                    end
                end
            end
        end
        a(TSend:365)=CurTS(Obs,6);
    else  % if it is not based on slow walking and slops
        tx1=1;
        y1=CurTS(1,6);
        for i=1:Obs;
            tx2=CurTS(i,4);
            y2=CurTS(i,6);
            a(tx1:tx2-1)=y1;
            tx1=tx2;
            y1=y2;
        end
      a(tx2:end)=y2;
    end
 
    
    %cond: time     : from TSstart=2 to TSend=340
    if strcmp(cond,'time')
        TS=TSstart;
        TE=TSend;
    end
    
    %cond: realSeq     : from TSstart=2 to x<SinceTime
    if strcmp(cond,'realSeq')
        x= CurTS(:,4)<=ToTime;
        x=CurTS(x,:);
        x=max(x(:,4));
        TS=TSstart;
        TE=x;
    end
    
    %cond: timeSeq  : from TSstart to ToTime
    if strcmp(cond,'timeSeq')
        TS=TSstart;
        TE=ToTime;
    end
    
    %cond: sequence : from SinceTime to ToTime
    if strcmp(cond,'sequence')
        TS=SinceTime;
        TE=ToTime;
        
    end
    
    a=a(TS:TE);
    if length(a)>10 && std(a) ~=0
        an =(a-mean(a))/std(a);
        inx=length(org_traj)+1;
        nor_traj{inx}=an;
        org_traj{inx}=a;
        t_traj{inx}=(TS:1:TE);
    end
end