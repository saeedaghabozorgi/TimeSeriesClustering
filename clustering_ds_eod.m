function clustering_ds_eod
% make_companies();
% make_train_ds();
% calculate_DTW_DIST();
load ('..\dataset eod\train data\KLSE_DATA.mat','train_data');
rows=size(train_data,1);
[nor_traj_raw,~]=z_normalize(1,rows,1,240,train_data);

calculate_DTW_matrix_paralele(nor_traj,'..\dataset eod\train data\KLSE_dismat_DTW');

dis_mtx=Mtx_Distance(nor_traj_raw,nor_traj_raw,'same','Norm','dis_method','Euclid');
k=round(sqrt(2*rows));
details=clustering_Hybrid_3Level(nor_traj_raw,k,randi(k,[rows,1]),[]);


end

function calculate_DTW_DIST()

end

function calculate_ED_DIST()
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
     make_figure(traj,file_name{k-2});
    rrr{k-2}= strrep(file_name{k-2},'.mat','')
end
save('EOD_DATA','train_data');
save('COMPANY_NAME','rrr')

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

function make_companies()
dd={'i-H.txt','i-I.txt','i-J.txt','i-K.txt','i-L.txt','i-M.txt','i-N.txt','i-O.txt','i-P.txt','i-Q.txt','i-R.txt','i-S.txt','i-T.txt','i-U.txt','i-V.txt','i-W.txt','i-X.txt','i-Y.txt','i-Z.txt'}
for i=1:length(dd)
    xx=dd{i}
    read_file_save(xx);
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


function calculate_DTW_matrix_paralele(nor_traj,path)
matlabpool open 8
poolSize = matlabpool('size');
if poolSize == 0
   error('parallel:demo:poolClosed', ...
        'This demo needs an open MATLAB pool to run.');
end
fprintf('This demo is running on %d MATLABPOOL workers.\n', ...
    matlabpool('size'));

x=nor_traj;
data_n = length(x);
dismat=zeros(data_n,data_n);
parfor  i = 1:data_n
    disp(num2str(i));
    dismat(i,:)=rowcalc(i,data_n,x);
end
dismat=dismat+dismat';
save(path, 'dismat')
matlabpool close
end

function dismatrix=rowcalc(row,data_n,x)
dismatrix=zeros(1,data_n);
    a=x{row};
    for j = row:data_n,
        if row ~= j
            b=x{j};
            dis=dis_dtw3(a,b,size(a,2));
            dismatrix(1, j) =dis;
        end
    end
end
