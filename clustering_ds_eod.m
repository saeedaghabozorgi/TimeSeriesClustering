function clustering_ds_eod

%read_file_save('i-D.txt');
make_train_ds()


end

function make_train_ds()
foldpath='C:\Users\Saeed\Documents\Google project\dataset eod\companies2010\';
nameFiles = dir(foldpath);
x=[];
for k=3:length(nameFiles)
    file_name{k-2} = nameFiles(k).name;
    file_path{k-2}=[foldpath,file_name{k-2}];
    ff=file_path{k-2};
    D=load(ff,'dd');
    traj=D.dd(1:240,4)';
    train_data(k-2,:)=traj;
%     make_figure(traj,file_name{k-2});
    rrr{k-2}= strrep(file_name{k-2},'.mat','')
end
rows=size(train_data,1);
[nor_traj_raw,~]=z_normalize(1,rows,1,240,train_data);


dis_mtx=Mtx_Distance(nor_traj_raw,nor_traj_raw,'same','Norm','dis_method','Euclid');
details=clustering_Hybrid_3Level(nor_traj_raw,16,randi(16,[rows,1]),[]);

end


function make_figure(traj,filename)
h = figure(1); 
plot( traj);
filename=['C:\Users\Saeed\Documents\Google project\dataset eod\plot2010\',filename,'.png'];
print(h,'-dpng',filename);
end

function [nor_traj,t_traj]=z_normalize(FromUser,ToUser,SinceTime,ToTime,data)
org_traj=[];  % orginal traj
nor_traj=[]; % normalized orginal traj
t_traj=[];

for z=FromUser:ToUser
    
    TS=SinceTime;
    TE=ToTime;
    a=data(z,TS:TE);
    an =(a-mean(a))/std(a);
    inx=length(org_traj)+1;
    nor_traj{inx}=an;
    org_traj{inx}=a;
    t_traj{inx}=(1:1:TE-TS+1);
end
end


function read_file_save(text_file)
%# Read the XLS file
folderpath='C:\Users\Saeed\Documents\Google project\dataset eod\';
read_data = importdata([folderpath,'raw data\',text_file]);
D=read_data.data;
txt=read_data.textdata(2:end,:);
%[D txt] = xlsread([folderpath,'A.csv']);
data=D(1,:);
dt=[];
j=1;
datem(j,1)=txt(1,1);
datem(j,2)=txt(1,3);
for i=2:length(txt)
    
    if(strcmp(txt(i,1),txt(i-1,1)))
        dt=[dt;datenum(txt(i,3))];
        data=[data;D(i,:)];
    else
        
        xx=find(dt>datenum('1/1/2010') & dt<datenum('1/1/2011'));
       % TS{j}=data(:,4)';
        datem(j,1)=txt(i-1,1);
        datem(j,3)=txt(i-1,3);
        filename=strcat(folderpath,'companies\',txt(i-1,1),'.mat');
        save(char(filename),'data');
        if ( length(xx)>200 & datenum(datem(j,2))> datenum('1/1/2011') & datenum(datem(j,3))< datenum('1/1/2010') )
            dd=data(xx,:);
            filename=strcat(folderpath,'companies2010\',txt(i-1,1),'.mat');
            save(char(filename),'dd');
        end
        data=D(i,:);
        dt=[];
        j=j+1
        datem(j,2)=txt(i,3);
    end
end
        filename=strcat(folderpath,'info\',text_file,'.mat');
        save(char(filename),'datem');
end